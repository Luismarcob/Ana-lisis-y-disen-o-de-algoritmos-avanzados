#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using namespace::std;

void fileRead(vector<string> &transmissionOne, vector<string> &transmissionTwo, string &codeOne, string &codeTwo, string &codeThree){
    ifstream transmissionOneFile("transmission1.txt");
    ifstream transmissionTwoFile("transmission2.txt");
    ifstream codeOneFile("mcode1.txt");
    ifstream codeTwoFile("mcode2.txt");
    ifstream codeThreeFile("mcode3.txt");

    string currentLine;

    while(getline(transmissionOneFile, currentLine)){
        transmissionOne.push_back(currentLine);
    }

    while(getline(transmissionTwoFile, currentLine)){
        transmissionTwo.push_back(currentLine);
    }

    while(getline(codeOneFile, currentLine)){
        codeOne = currentLine;
    }

    while(getline(codeTwoFile, currentLine)){
        codeTwo = currentLine;
    }

    while(getline(codeThreeFile, currentLine)){
        codeThree = currentLine;
    }

    transmissionOneFile.close();
    transmissionTwoFile.close();
    codeOneFile.close();
    codeTwoFile.close();
    codeThreeFile.close();
}

void printFiles(vector<string> transmissionOne, vector<string> transmissionTwo, string codeOne, string codeTwo, string codeThree){
    cout << "TRANSMISSION 1" << endl;

    for(int i = 0; i < transmissionOne.size(); i++){
        cout << transmissionOne[i] << endl;
    }

    cout << "TRANSMISSION 2" << endl;
    for(int i = 0; i < transmissionTwo.size(); i++){
        cout << transmissionTwo[i] << endl;
    }

    cout << "CODE 1: " << codeOne << endl;
    cout << "CODE 2: " << codeTwo << endl;
    cout << "CODE 3: " << codeTwo << endl;

}


vector<int> preKMP(string code){
    int m = code.length(); //Obtenemos el largo de nuestro string
    vector<int> lps(m,0); //Vector LPS inicializado en 0

    int j = 0;
    int i = 1;

    while(i < m){
        if(code[i] == code[j]){
            lps[i] = j+1;
            i++;
            j++;
        } 
        else{ 
            if (j > 0) {
                j = lps[j-1];
            } 
            else {
                lps[i] = 0;
                i++;
            }
        }
    }
    return lps;
}

bool KMP(vector<string> transmission, string code) {
    int n = transmission.size();
    int m = code.length();

    // Generar el vector LPS
    vector<int> lps = preKMP(code);

    // Para cada línea de la transmisión, buscamos el patrón
    for (int t = 0; t < n; t++) {
        string text = transmission[t];
        int i = 0; // índice para text
        int j = 0; // índice para code

        while (i < text.length()) {
            if (code[j] == text[i]) {
                i++;
                j++;
            }

            if (j == m) {
                // Se encontró el patrón en la posición (i-j)
                return true;
                j = lps[j - 1]; // Continuar buscando otras ocurrencias
            } else if (i < text.length() && code[j] != text[i]) {
                // Desajuste después de j coincidencias
                if (j != 0) {
                    j = lps[j - 1];
                } else {
                    i++;
                }
            }
        }
    }
    return false;
}


string preManacher(string code){
    string newString = "";

    for(int i = 0; i < code.size(); i++){
        newString.append("#");
        newString.append(string(1,code[i]));
    } 
    newString.append("#");
    return newString;
}

void manacher(vector<string> transmission) {

    for (int t = 0; t < transmission.size(); t++) {
        // Procesar el input
        string T = preManacher(transmission[t]);
        int n = T.length();
        vector<int> P(n, 0); // Arreglo para el largo del palindromo en cada posicion
        int C = 0, R = 0; // Border izquierdo y derecho

        for (int i = 1; i < n - 1; i++) {
            int mirror = 2 * C - i; // Posicion espejo

            // Si i estan en el  borde derecho
            if (i < R) {
                P[i] = min(R - i, P[mirror]);
            }

            // Intentar expandir el palindromo
            while (T[i + P[i] + 1] == T[i - P[i] - 1]) {
                P[i]++;
            }

            // Actualizar el palindromo
            if (i + P[i] > R) {
                C = i;
                R = i + P[i];
            }
        }

        // Encontrar el palindromo mas largo
        int maxLen = 0;
        int centerIndex = 0;
        for (int i = 1; i < n - 1; i++) {
            if (P[i] > maxLen) {
                maxLen = P[i];
                centerIndex = i;
            }
        }

        cout << "\t" << "Linea: " << t << " Inicio: " << (centerIndex - maxLen)/2 << " Fin: " << ((centerIndex - maxLen)/2) + maxLen - 1 << endl;


    }
}

void LCS(vector<string> transmissionOne, vector<string> transmissionTwo) { //Longest Common Substring
    string transmissionOneString = "";
    string transmissionTwoString = "";

    // Concatena las transmisiones
    for(int i = 0; i < transmissionOne.size(); i++) {
        transmissionOneString.append(transmissionOne[i]);
    }

    for(int i = 0; i < transmissionTwo.size(); i++) {
        transmissionTwoString.append(transmissionTwo[i]);
    }

    int n = transmissionOneString.length();
    int m = transmissionTwoString.length();
    int maxSize = 0; // Tamaño del substring mas largo
    int endIndex = 0; // Indice final del substring mas largo

    // Creacion de la table donde guardaremos el largo de las substrings
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));

    // Llenar la tabla
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            if (transmissionOneString[i - 1] == transmissionTwoString[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1] + 1;
                if (dp[i][j] > maxSize) {
                    maxSize = dp[i][j];
                    endIndex = i - 1;
                }
            } else {
                dp[i][j] = 0;
            }
        }
    }


    // Calculamos las posiciones
    int startIndex = endIndex - maxSize + 2; 
    int endIndexFinal = endIndex + 1; 

    cout << "\t" << startIndex << " " << endIndexFinal << endl;

    // Longest common substring
    string longestCommonSubstring = transmissionOneString.substr(endIndex - maxSize + 1, maxSize);

    cout << "\t" << "Longest Common Substring: " << longestCommonSubstring << endl;
}


int main(){
    vector<string> transmissionOne;
    vector<string> transmissionTwo;
    string codeOne;
    string codeTwo;
    string codeThree;

    fileRead(transmissionOne,transmissionTwo,codeOne,codeTwo,codeThree);
    vector<int> codeOneLps = preKMP(codeOne);

    //Buscar el codigo malicioso - Transmision 1
    cout << "Buscando Codigo Malicioso 1 en transmision 1" << endl;
    cout << "\t" << boolalpha << KMP(transmissionOne, codeOne) << endl; 

    cout << "Buscando Codigo Malicioso 2 en transmision 1" << endl;
    cout << "\t" << KMP(transmissionOne, codeTwo) << endl; 

    cout << "Buscando Codigo Malicioso 3 en transmision 1" << endl;
    cout << "\t" << KMP(transmissionOne, codeThree) << endl; 

    //Buscar el codigo malicioso - Transmision 2
    cout << "Buscando Codigo Malicioso 1 en transmision 2" << endl;
    cout << "\t" << KMP(transmissionTwo, codeOne) << endl; 

    cout << "Buscando Codigo Malicioso 2 en transmision 2" << endl;
    cout << "\t" << KMP(transmissionTwo, codeTwo) << endl; 

    cout << "Buscando Codigo Malicioso 3 en transmision 2" << endl;
    cout << "\t" << KMP(transmissionTwo, codeThree) << endl; 


    //Palindromo mas largo;

    cout << "Palindromo mas largo transmision 1" << endl;
    manacher(transmissionOne);

    cout << "Palindromo mas largo transmision 2" << endl;
    manacher(transmissionTwo);

    cout << "Buscando substring mas largo" << endl;
    LCS(transmissionOne, transmissionTwo);

    //cout << preManacher(codeOne);

    //Funcion para imprimir y aseguranos de la integridad de los inputs
    //printFiles(transmissionOne,transmissionTwo,codeOne,codeTwo,codeThree);

    /*
    for(int i = 0; i < codeOneLps.size(); i++){
        cout << codeOneLps[i] << endl;
    }
    */



    return 0;
}