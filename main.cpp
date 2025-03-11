#include <iostream>
#include <vector>
#include <limits.h>
#include <algorithm>
#include <queue>
#include <cmath>
#include <set>
#include "voronoi.h"
using namespace std;

const int VAL_INFINITO = INT_MAX; //Valor "infinito"

// Algoritmo de Prim. Complejidad: O(V^2).
vector<pair<int, int>> primMST(const vector<vector<int>>& GRAFO, int N) {
    vector<int> pesoMin(N, VAL_INFINITO); //Vector de pesos minimos        
    vector<bool> enMST(N, false);    //Vector de nodos visitados
    vector<int> nodoAnterior(N, -1);     //Vector de padres
    
    //Nodo raiz
    pesoMin[0] = 0;               
    nodoAnterior[0] = -1;              

    for (int count = 0; count < N - 1; ++count) { //Recorremos todos los nodos. O(V)
        int u = -1;

        // Encontramos el nodo no visitado con el peso menor
        for (int i = 0; i < N; ++i) {
            if (!enMST[i] && (u == -1 || pesoMin[i] < pesoMin[u])) { //Si no está en el MST y el peso es menor o es el primer nodo
                u = i; //Actualizamos el nodo
            }
        }

        enMST[u] = true; //Visitado

        // Actualizamos el peso de los nodos adjacentes
        for (int v = 0; v < N; ++v) {
            if (GRAFO[u][v] && !enMST[v] && GRAFO[u][v] < pesoMin[v]) { //Si no está en el MST y el peso es menor
                pesoMin[v] = GRAFO[u][v];
                nodoAnterior[v] = u;
            }
        }
    }

    //Lista de arcos para el MST
    vector<pair<int, int>> rFinalMST;
    for (int i = 1; i < N; ++i) {
        rFinalMST.push_back({nodoAnterior[i], i});
    }
    return rFinalMST;
}

// Algoritmo Held-Karp para el TSP usando programación dinámica. Complejiad: O(2^N * N^2).
int tsp(int posNodo, int BitMaskVisitado, int N, const vector<vector<int>>& GRAFO, vector<vector<int>>& dpTabla) {
    if (BitMaskVisitado == (1 << N) - 1) {
        return GRAFO[posNodo][0]; 
    }
    if (dpTabla[posNodo][BitMaskVisitado] != -1) {
        return dpTabla[posNodo][BitMaskVisitado];
    }
    
    int costoMin = VAL_INFINITO;
    for (int city = 0; city < N; ++city) {
        if ((BitMaskVisitado & (1 << city)) == 0) { 
            int nuevoCosto = GRAFO[posNodo][city] + tsp(city, BitMaskVisitado | (1 << city), N, GRAFO, dpTabla);
            costoMin = min(costoMin, nuevoCosto);
        }
    }
    return dpTabla[posNodo][BitMaskVisitado] = costoMin;
}

//Reconstruir la ruta óptima
void reconstructRoute(int posNodo, int BitMaskVisitado, int N, const vector<vector<int>>& GRAFO, vector<vector<int>>& dpTabla, vector<int>& path) {
    if (BitMaskVisitado == (1 << N) - 1) {
        return;
    }

    int costoMin = VAL_INFINITO;
    int nextCity = -1;

    for (int city = 0; city < N; ++city) {
        if ((BitMaskVisitado & (1 << city)) == 0) {  
            int nuevoCosto = GRAFO[posNodo][city] + dpTabla[city][BitMaskVisitado | (1 << city)];
            if (nuevoCosto < costoMin) {
                costoMin = nuevoCosto;
                nextCity = city;
            }
        }
    }
    
    path.push_back(nextCity);
    reconstructRoute(nextCity, BitMaskVisitado | (1 << nextCity), N, GRAFO, dpTabla, path);
}

//Búsqueda en anchura para encontrar un camino de aumento. Complejidad: O(V^2).
bool bfs(const vector<vector<int>>& capacidad, vector<vector<int>>& flow, vector<int>& nodoAnterior, int N, int source, int sink) {
    vector<bool> BitMaskVisitado(N, false);
    queue<int> q;
    q.push(source);
    BitMaskVisitado[source] = true;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int v = 0; v < N; ++v) {
            if (!BitMaskVisitado[v] && capacidad[u][v] - flow[u][v] > 0) {
                q.push(v);
                BitMaskVisitado[v] = true;
                nodoAnterior[v] = u;
                if (v == sink) return true;
            }
        }
    }
    return false;
}

// Algoritmo de Ford-Fulkerson para calcular el flujo máximo. Complejidad: O(V * E^2).
int fordFulkerson(const vector<vector<int>>& capacidad, int N, int source, int sink) {
    vector<vector<int>> flow(N, vector<int>(N, 0)); 
    vector<int> nodoAnterior(N);

    int maxFlujo = 0;

    // Mientras haya un camino de aumento
    while (bfs(capacidad, flow, nodoAnterior, N, source, sink)) {
        // Encuentra el flujo máximo posible en el camino de aumento
        int pathFlow = VAL_INFINITO; 
        for (int v = sink; v != source; v = nodoAnterior[v]) { //Recorremos el camino de aumento
            int u = nodoAnterior[v]; //Nodo anterior
            pathFlow = min(pathFlow, capacidad[u][v] - flow[u][v]); //Flujo máximo
        }

        // Actualiza los flujos a lo largo del camino de aumento
        for (int v = sink; v != source; v = nodoAnterior[v]) {
            int u = nodoAnterior[v];
            flow[u][v] += pathFlow;
            flow[v][u] -= pathFlow; 
        }

        maxFlujo += pathFlow;
    }
    return maxFlujo;
}


//La siguente seccion ha sido obtenida de https://www.cs.hmc.edu/~mbrubeck/voronoi.html
//Creditos a Matt Brubeck

#pragma region voronoi //Complejidad: O(N log N)

priority_queue<point,  vector<point>,  gt> points; // site events
priority_queue<event*, vector<event*>, gt> events; // circle events

void process_point()
{
   // Get the next point from the queue.
   point p = points.top();
   points.pop();

   // Add a new arc to the parabolic front.
   front_insert(p);
}

void process_event()
{
   // Get the next event from the queue.
   event *e = events.top();
   events.pop();

   if (e->valid) {
      // Start a new edge.
      seg *s = new seg(e->p);

      // Remove the associated arc from the front.
      arc *a = e->a;
      if (a->prev) {
         a->prev->next = a->next;
         a->prev->s1 = s;
      }
      if (a->next) {
         a->next->prev = a->prev;
         a->next->s0 = s;
      }

      // Finish the edges before and after a.
      if (a->s0) a->s0->finish(e->p);
      if (a->s1) a->s1->finish(e->p);

      // Recheck circle events on either side of p:
      if (a->prev) check_circle_event(a->prev, e->x);
      if (a->next) check_circle_event(a->next, e->x);
   }
   delete e;
}

void front_insert(point p)
{
   if (!root) {
      root = new arc(p);
      return;
   }

   // Find the current arc(s) at height p.y (if there are any).
   for (arc *i = root; i; i = i->next) {
      point z, zz;
      if (intersect(p,i,&z)) {
         // New parabola intersects arc i.  If necessary, duplicate i.
         if (i->next && !intersect(p,i->next, &zz)) {
            i->next->prev = new arc(i->p,i,i->next);
            i->next = i->next->prev;
         }
         else i->next = new arc(i->p,i);
         i->next->s1 = i->s1;

         // Add p between i and i->next.
         i->next->prev = new arc(p,i,i->next);
         i->next = i->next->prev;

         i = i->next; // Now i points to the new arc.

         // Add new half-edges connected to i's endpoints.
         i->prev->s1 = i->s0 = new seg(z);
         i->next->s0 = i->s1 = new seg(z);

         // Check for new circle events around the new arc:
         check_circle_event(i, p.x);
         check_circle_event(i->prev, p.x);
         check_circle_event(i->next, p.x);

         return;
      }
   }

   // Special case: If p never intersects an arc, append it to the list.
   arc *i;
   for (i = root; i->next; i=i->next) ; // Find the last node.

   i->next = new arc(p,i);
   // Insert segment between p and i
   point start;
   start.x = X0;
   start.y = (i->next->p.y + i->p.y) / 2;
   i->s1 = i->next->s0 = new seg(start);
}

// Look for a new circle event for arc i.
void check_circle_event(arc *i, double x0)
{
   // Invalidate any old event.
   if (i->e && i->e->x != x0)
      i->e->valid = false;
   i->e = NULL;

   if (!i->prev || !i->next)
      return;

   double x;
   point o;

   if (circle(i->prev->p, i->p, i->next->p, &x,&o) && x > x0) {
      // Create new event.
      i->e = new event(x, o, i);
      events.push(i->e);
   }
}

// Find the rightmost point on the circle through a,b,c.
bool circle(point a, point b, point c, double *x, point *o)
{
   // Check that bc is a "right turn" from ab.
   if ((b.x-a.x)*(c.y-a.y) - (c.x-a.x)*(b.y-a.y) > 0)
      return false;

   // Algorithm from O'Rourke 2ed p. 189.
   double A = b.x - a.x,  B = b.y - a.y,
          C = c.x - a.x,  D = c.y - a.y,
          E = A*(a.x+b.x) + B*(a.y+b.y),
          F = C*(a.x+c.x) + D*(a.y+c.y),
          G = 2*(A*(c.y-b.y) - B*(c.x-b.x));

   if (G == 0) return false;  // Points are co-linear.

   // Point o is the center of the circle.
   o->x = (D*E-B*F)/G;
   o->y = (A*F-C*E)/G;

   // o.x plus radius equals max x coordinate.
   *x = o->x + sqrt( pow(a.x - o->x, 2) + pow(a.y - o->y, 2) );
   return true;
}

// Will a new parabola at point p intersect with arc i?
bool intersect(point p, arc *i, point *res)
{
   if (i->p.x == p.x) return false;

   double a,b;
   if (i->prev) // Get the intersection of i->prev, i.
      a = intersection(i->prev->p, i->p, p.x).y;
   if (i->next) // Get the intersection of i->next, i.
      b = intersection(i->p, i->next->p, p.x).y;

   if ((!i->prev || a <= p.y) && (!i->next || p.y <= b)) {
      res->y = p.y;

      // Plug it back into the parabola equation.
      res->x = (i->p.x*i->p.x + (i->p.y-res->y)*(i->p.y-res->y) - p.x*p.x)
                / (2*i->p.x - 2*p.x);

      return true;
   }
   return false;
}

// Where do two parabolas intersect?
point intersection(point p0, point p1, double l)
{
   point res, p = p0;

   if (p0.x == p1.x)
      res.y = (p0.y + p1.y) / 2;
   else if (p1.x == l)
      res.y = p1.y;
   else if (p0.x == l) {
      res.y = p0.y;
      p = p1;
   } else {
      // Use the quadratic formula.
      double z0 = 2*(p0.x - l);
      double z1 = 2*(p1.x - l);

      double a = 1/z0 - 1/z1;
      double b = -2*(p0.y/z0 - p1.y/z1);
      double c = (p0.y*p0.y + p0.x*p0.x - l*l)/z0
               - (p1.y*p1.y + p1.x*p1.x - l*l)/z1;

      res.y = ( -b - sqrt(b*b - 4*a*c) ) / (2*a);
   }
   // Plug back into one of the parabola equations.
   res.x = (p.x*p.x + (p.y-res.y)*(p.y-res.y) - l*l)/(2*p.x-2*l);
   return res;
}

void finish_edges()
{
   // Advance the sweep line so no parabolas can cross the bounding box.
   double l = X1 + (X1-X0) + (Y1-Y0);

   // Extend each remaining segment to the new parabola intersections.
   for (arc *i = root; i->next; i = i->next)
      if (i->s1)
         i->s1->finish(intersection(i->p, i->next->p, l*2));
}

void print_output()
{
   // Bounding box coordinates.
   cout << "Coordenadas de la caja delimitadora X0, X1, Y0, Y1:\n";
   cout << "\t" << X0 << " "<< X1 << " " << Y0 << " " << Y1 << endl;

   // Each output segment in four-column format.
   cout << "Segmentos de salida: X0, X1, Y0, Y1\n";
   vector<seg*>::iterator i;
   for (i = output.begin(); i != output.end(); i++) {
        point p0 = (*i)->start;
        point p1 = (*i)->end;
        cout << "\t" << p0.x << " " << p0.y << " " << p1.x << " " << p1.y << endl;
   }
}
#pragma endregion


int main() {
    cout << "Introduzca el numero de ciudades: ";
    int N;
    cin >> N;

    // Leer la matriz de distancias
    vector<vector<int>> GRAFO(N, vector<int>(N));
    for (int i = 0; i < N; ++i) {
        cout << "Introduzca las distancias de la ciudad " << char('A' + i) << " a las demas ciudades." << endl;
        for (int j = 0; j < N; ++j) {
            cin >> GRAFO[i][j];
        }
    }

    // Leer la matriz de capacidades de flujo
    vector<vector<int>> capacidad(N, vector<int>(N));
    for (int i = 0; i < N; ++i) {
        cout << "Introduzca las capacidades de flujo de la ciudad " << char('A' + i) << " a las demas ciudades." << endl;
        for (int j = 0; j < N; ++j) {
            cin >> capacidad[i][j];
        }
    }

    
    // Leer las coordenadas de las centrales
    // Input de voronoi obtenido de: https://www.cs.hmc.edu/~mbrubeck/voronoi.html
    //Creditos a Matt Brubeck
    point p;
    char ignore;
    for(int i = 0; i < N; ++i) {
        cout << "Introduzca las coordenadas de la central " << char('A' + i) << ": ";
        cin >> ignore >> p.x >> ignore >> p.y >> ignore; // Read format (x,y)
        points.push(p);

        // Mantener el tamaño de la caja delimitadora
        if (p.x < X0) X0 = p.x;
        if (p.y < Y0) Y0 = p.y;
        if (p.x > X1) X1 = p.x;
        if (p.y > Y1) Y1 = p.y;
    }
    // Agregar márgenes a la caja delimitadora.
    double dx = (X1-X0+1)/5.0, dy = (Y1-Y0+1)/5.0;
    X0 -= dx;  X1 += dx;  Y0 -= dy;  Y1 += dy;

    // Procesar las colas; seleccionar el elemento superior con la coordenada x más pequeña.
    while (!points.empty())
        if (!events.empty() && events.top()->x <= points.top().x)
            process_event();
        else
            process_point();

    // Después de procesar todos los puntos, procesar los eventos de círculo restantes.
    while (!events.empty())
        process_event();

    //MST usando el algoritmo de Prim
    vector<pair<int, int>> rFinalMST = primMST(GRAFO, N);

    //lista de arcos que forman el MST

    cout << "" << endl;
    cout << "Cableado optimo entre colonias:" << endl;
    for (const auto& edge : rFinalMST) {
        char u = 'A' + edge.first;
        char v = 'A' + edge.second;
        cout << "(" << u << ", " << v << ")\n";
    }

    //ruta óptima usando TSP
    vector<vector<int>> dpTabla(N, vector<int>(1 << N, -1));
    int minPathCost = tsp(0, 1, N, GRAFO, dpTabla);

    // Reconstruimos la ruta del TSP
    vector<int> path = {0}; 
    reconstructRoute(0, 1, N, GRAFO, dpTabla, path);

    // Mostramos la ruta óptima
    cout << "\nRuta optima para distribucion:\n";
    for (int i = 0; i < path.size(); ++i) {
        cout << char('A' + path[i]) << " -> ";
    }
    cout << "A" << endl; 
    cout << "Costo minimo de la ruta: " << minPathCost << endl;

    // Calculamos el flujo máximo entre la colonia 0 (inicio) y la colonia N-1 (final)
    int maxFlujo = fordFulkerson(capacidad, N, 0, N - 1);

    cout << "\nEl flujo maximo de informacion entre la colonia A (origen) y la colonia " 
         << char('A' + N - 1) << " (destino) es: " << maxFlujo << "\n" << endl;

         
    finish_edges(); // Clean up dangling edges.
    print_output(); // Output the voronoi diagram.

    return 0;
}
