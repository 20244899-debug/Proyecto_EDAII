#include <iostream>
#include <fstream> //esto es para poder abrir archivo
#include <string>
#include <vector>
#include <cmath> //Para el exponente en la parte d)
#include <iomanip> //PREGUNTAR PARA QUE SIRVE


using namespace std;//esta es para no tener que escribir std::, lo de la mañana
//esta es la parte para que se pueda ingresar el archivo
struct DNASequence {
    string id;
    string sequence;
};

bool leerFasta(string ruta, vector<DNASequence>& lista/*aqui lo va a guardar en una lista*/) {
    ifstream archivo(ruta.c_str());//el doc que se ingreso se va a llamar archivo
    if (!archivo.is_open()) return false;

    string linea;//guarda temporalmente una linea
    DNASequence actual;
    bool hayActual = false;

    while (getline(archivo, linea)) {
        if (!linea.empty() && linea[linea.size() - 1] == '\r') linea.erase(linea.size() - 1);
        if (linea.empty()) continue;//esto es para ignorar lineas vacia

        if (linea[0] == '>') {//es para ver si es el ID de la secuencia
            if (hayActual) lista.push_back(actual);//es para cuando se inicie una nueva secuencia se agregue al actual
            actual.id = linea.substr(1);//es para quitar > del nombre del id
            actual.sequence = "";
            hayActual = true;
        } else {
            actual.sequence += linea;//esto es para la secuencia, si esta en otra linea unirlas
        }
    }
    if (hayActual) lista.push_back(actual);//esto es para guardar la ultima secuencia
    return true;
}

//esta es una funcion para sacar si es primo (para despues)
bool esPrimo(int n) {
    if (n < 2) return false;
    for (int i = 2; i * i <= n; i++)
        if (n % i == 0) return false;
    return true;
}
//funcion para buscar el siguiente primo
int siguientePrimo(int n) {
    while (!esPrimo(n)) n++;
    return n;
}

//esta es la funcion hash
int HashFunction(string key, int a, int length) {
    unsigned long long h = 0;
    for (int i = 0; i < (int)key.size(); i++)//recorre toda la secuencia de letras
        h = (h * a + (unsigned char)key[i]) % length;
    return (int)h;
}

// ---------- PARTE B, TABLA 1: LINEAR PROBING ----------
struct Bin {
    DNASequence dato;
    bool ocupado;
};

struct HashTable {
    vector<Bin> bins;//este es un vector que guarda las secuencias
    int length;//cuantos espacios existen
    int num_keys;//cuantos espacios estan ocupados
};

void crearTabla(HashTable& ht, int tam) {
    ht.length = tam;
    ht.num_keys = 0;
    ht.bins.assign(tam, Bin());//con el tamaño crea la cant de cajones
    for (int i = 0; i < tam; i++) ht.bins[i].ocupado = false;//pone que cajones estan vacios
}

int HashTableBusqueda(HashTable& ht, string key) {
    int index = HashFunction(key, 37, ht.length);//que posicion de la tabla buscar
    int conteo = 0;
    while (ht.bins[index].ocupado && conteo != ht.length) {
        if (ht.bins[index].dato.id == key) return index;
        index = (index + 1) % ht.length;
        conteo = conteo + 1;
    }
    return -1;
}
//esta funcion es para aumentar el tamaño de la tabla
void redimensionar(HashTable& ht) {
    vector<Bin> vieja = ht.bins;
    crearTabla(ht, siguientePrimo(2 * (ht.num_keys + 1)));

    for (int i = 0; i < (int)vieja.size(); i++) {
        if (vieja[i].ocupado) {
            int index = HashFunction(vieja[i].dato.id, 37, ht.length);
            while (ht.bins[index].ocupado) index = (index + 1) % ht.length;
            ht.bins[index] = vieja[i];
            ht.num_keys = ht.num_keys + 1;
        }
    }
}   // <-- esta llave cierra la funcion

bool HashTableInsertar(HashTable& ht, DNASequence s) {
    int pos = HashTableBusqueda(ht, s.id);
    if (pos != -1) {
        ht.bins[pos].dato = s;
        return true;
    }
    if ((double)(ht.num_keys + 1) / ht.length > 0.7) redimensionar(ht);

    int index = HashFunction(s.id, 37, ht.length);
    int conteo = 0;
    while (ht.bins[index].ocupado && conteo != ht.length) {
        index = (index + 1) % ht.length;
        conteo = conteo + 1;
    }
    if (conteo == ht.length) return false;

    ht.bins[index].dato = s;
    ht.bins[index].ocupado = true;
    ht.num_keys = ht.num_keys + 1;
    return true;
}

// ---------- PARTE B, TABLA 2: DOUBLE HASHING (k-mers #el terrorrr) ----------
struct BinKmer {
    string key;//es el kmer como ATG
    int valor;//es cuantas veces aparece eso
    bool ocupado;
};

struct HashTableKmer {
    vector<BinKmer> bins;
    int length;      // 4^k
    int num_keys;
    int k;
};

void crearTablaKmer(HashTableKmer& ht, int k) {
    ht.k = k;
    ht.length = 1;
    for (int i = 0; i < k; i++) ht.length = ht.length * 4;
    ht.num_keys = 0;
    ht.bins.assign(ht.length, BinKmer());
    for (int i = 0; i < ht.length; i++) {
        ht.bins[i].valor = 0;
        ht.bins[i].ocupado = false;
    }
}

// h1: indice inicial.  h2: salto (impar, porque length es potencia de 2)
int hash1(string key, int length) { return HashFunction(key, 37, length); }
int hash2(string key, int length) { return 2 * HashFunction(key, 41, length / 2) + 1; }

bool KmerInsertar(HashTableKmer& ht, string key) {
    int index = hash1(key, ht.length);
    int salto = hash2(key, ht.length);
    int conteo = 0;

    while (conteo != ht.length) {
        if (!ht.bins[index].ocupado) {
            ht.bins[index].key = key;
            ht.bins[index].valor = 1;
            ht.bins[index].ocupado = true;
            ht.num_keys = ht.num_keys + 1;
            return true;
        }
        if (ht.bins[index].key == key) {
            ht.bins[index].valor = ht.bins[index].valor + 1;
            return true;
        }
        index = (index + salto) % ht.length;      // double hashing
        conteo = conteo + 1;
    }
    return false;
}

int KmerBusqueda(HashTableKmer& ht, string key) {
    int index = hash1(key, ht.length);
    int salto = hash2(key, ht.length);
    int conteo = 0;

    while (ht.bins[index].ocupado && conteo != ht.length) {
        if (ht.bins[index].key == key) return ht.bins[index].valor;
        index = (index + salto) % ht.length;
        conteo = conteo + 1;
    }
    return 0;
}

// Ventana de k letras avanzando de 1 en 1; se ignoran k-mers con letras distintas de ACGT
void construirKmers(HashTableKmer& ht, string seq) {
    for (int i = 0; i + ht.k <= (int)seq.size(); i++) {
        string kmer = seq.substr(i, ht.k);
        bool valido = true;
        for (int j = 0; j < (int)kmer.size(); j++)
            if (kmer[j] != 'A' && kmer[j] != 'C' && kmer[j] != 'G' && kmer[j] != 'T') valido = false;
        if (valido) KmerInsertar(ht, kmer);
    }
}

// ---------- PARTE C: Max-Heap ----------

struct KmerCant {
    string kmer;
    int conteo; //guarda la cantidad de veces que aparece
};

//Posiciones en el heap
int padre(int i)  { return (i - 1) / 2; }
int hijoIzq(int i) { return 2 * i + 1; }
int hijoDer(int i) { return 2 * i + 2; }


//Solo se hace uso del heapifyDown
void heapifyDown(vector<KmerCant>& heap, int tamano, int i){
    int mayor = i;
    int L = hijoIzq(i);
    int R = hijoDer(i);

    if (L < tamano && heap[L].conteo > heap[mayor].conteo) mayor = L;
    if (R < tamano && heap[R].conteo > heap[mayor].conteo) mayor = R;

    if (mayor != i){
        swap(heap[i], heap[mayor]);
        heapifyDown(heap, tamano, mayor);
    }
}

//Las dos funciones se va a necesitar para hallar los k-mers más frecuentes

void construirHeap(vector<KmerCant>& heap) {
    for (int i = padre(heap.size()- 1); i >= 0; i--)
        heapifyDown(heap, heap.size(), i);
}

KmerCant extraerMax(vector<KmerCant>& heap) {
    KmerCant mayor = heap[0];
    heap[0] = heap.back();
    heap.pop_back();
    if (!heap.empty()) {
        heapifyDown(heap, heap.size(), 0);
    }
    return mayor;
}

void topNKmers(HashTableKmer& tabla, int n) {
    vector<KmerCant> heap;
    for (int i = 0; i < tabla.length; i++) {
        if (tabla.bins[i].ocupado) {
            KmerCant kc;
            kc.kmer = tabla.bins[i].key;
            kc.conteo = tabla.bins[i].valor;
            heap.push_back(kc);
        }
    }

    construirHeap(heap);

    cout << "Top " << n << " k-mers mas frecuentes:" << endl;
    for (int i = 0; i < n && !heap.empty(); i++) {
        KmerCant kc = extraerMax(heap);
        cout << i + 1 << ". " << kc.kmer << " -> " << kc.conteo << endl;
    }
}

// ---------- PARTE D: Grafico de secuencias de ADN (Matriz de adyacencia) ----------

//Para crearr una matriz necesitamos indices numéricos
int indiceBase(char c) {
    if (c == 'A') return 0;
    if (c == 'C') return 1;
    if (c == 'G') return 2;
    if (c == 'T') return 3;
    return -1; //Si es que coloca otra letra es inválida
}

//Sirve para convertir indices numericos a 2-mer
string nombre2mer(int c) {
    const char letras[4] = {'A', 'C', 'G', 'T'};
    string s = "";
    s += letras[c / 4];
    s += letras[c % 4];
    return s;
}

void matrizAdyacencia(string secuencia, double beta) {
    int m = secuencia.size() / 2; //Si es que es una longitud impar, se descarta la ultima letra
    vector<int> par(m);
    bool usado[16] = {false}; // 16 porque es las combinaciones posible con 2-mer


    for (int p = 0; p < m; p++) {
        int a = indiceBase(secuencia[2 * p]);
        int b = indiceBase(secuencia[2 * p + 1]);
        if (a == -1 || b == -1) par[p] = -1;
        else {
            par[p] = a * 4 + b;
            usado[par[p]] = true;
        }
    }

    double A[16][16] = {0};
    for (int i = 0; i < m; i++)
        for (int j = i + 1; j < m; j++)
            if (par[i] != -1 && par[j] != -1)
                A[par[i]][par[j]] += 1.0 / pow(j - i, beta);

    cout << fixed << setprecision(3);
    cout << setw(6) << " ";
    for (int c = 0; c < 16; c++)
        if (usado[c]) cout << setw(8) << nombre2mer(c);
    cout << endl;

    for (int f = 0; f < 16; f++) {
        if (!usado[f]) continue;
        cout << setw(6) << nombre2mer(f);
        for (int c = 0; c < 16; c++)
            if (usado[c]) cout << setw(8) << A[f][c];
        cout << endl;
    }
}

// ---------- MENU ----------
int main() {
    HashTable ht;
    crearTabla(ht, 10);
    HashTableKmer hk;
    bool hayKmers = false;
    string idActual = "";
    int opcion = -1;

    while (opcion != 0) {
        cout << "\n1. Cargar FASTA"
             << "\n2. Buscar secuencia por ID"
             << "\n3. Estadisticas tabla 1"
             << "\n4. Construir tabla de k-mers"
             << "\n5. Buscar un k-mer"
             << "\n6. Mostrar k-mers"
             << "\n7. Top n k-mers mas frecuentes"
             << "\n8. Matriz de adyacencia de una secuencia"
             << "\n0. Salir\nOpcion: ";
        cin >> opcion;
        if (cin.fail()) { cin.clear(); cin.ignore(1000, '\n'); opcion = -1; continue; }

        if (opcion == 1) {
            string ruta;
            cout << "Ruta del archivo: ";
            cin.ignore(1000, '\n');
            getline(cin, ruta);
            if (!ruta.empty() && ruta[0] == '"') ruta.erase(0, 1);
            if (!ruta.empty() && ruta[ruta.size() - 1] == '"') ruta.erase(ruta.size() - 1);

            vector<DNASequence> lista;
            if (leerFasta(ruta, lista)) {
                for (int i = 0; i < (int)lista.size(); i++) HashTableInsertar(ht, lista[i]);
                cout << "Se leyeron " << lista.size() << " secuencias." << endl;
            } else {
                cout << "No se pudo abrir el archivo: " << ruta << endl;
            }
        } else if (opcion == 2) {
            string id;
            cout << "ID: ";
            cin >> id;
            int pos = HashTableBusqueda(ht, id);
            if (pos == -1) cout << "No encontrada." << endl;
            else cout << ht.bins[pos].dato.id << " (" << ht.bins[pos].dato.sequence.size()
                      << " letras): " << ht.bins[pos].dato.sequence.substr(0, 60) << "..." << endl;
        } else if (opcion == 3) {
            cout << "Keys: " << ht.num_keys << ", bins: " << ht.length
                 << ", load factor: " << (double)ht.num_keys / ht.length << endl;
        } else if (opcion == 4) {
            string id;
            int k;
            cout << "ID de la secuencia: ";
            cin >> id;
            int pos = HashTableBusqueda(ht, id);
            if (pos == -1) { cout << "No encontrada." << endl; continue; }
            cout << "Valor de k (1 a 10): ";
            cin >> k;
            if (k < 1 || k > 10) { cout << "k invalido." << endl; continue; }
            crearTablaKmer(hk, k);
            construirKmers(hk, ht.bins[pos].dato.sequence);
            hayKmers = true;
            idActual = id;
            cout << hk.length << " bins, " << hk.num_keys << " k-mers distintos." << endl;
        } else if (opcion == 5) {
            if (!hayKmers) { cout << "Primero use la opcion 4." << endl; continue; }
            string kmer;
            cout << "K-mer de " << hk.k << " letras: ";
            cin >> kmer;
            cout << kmer << " aparece " << KmerBusqueda(hk, kmer) << " veces en " << idActual << endl;
        } else if (opcion == 6) {
            if (!hayKmers) { cout << "Primero use la opcion 4." << endl; continue; }
            int mostrados = 0;
            for (int i = 0; i < hk.length && mostrados < 20; i++)
                if (hk.bins[i].ocupado) {
                    cout << hk.bins[i].key << "\t" << hk.bins[i].valor << endl;
                    mostrados++;
                }
        } else if (opcion == 7) {
            if (!hayKmers) { cout << "Primero use la opcion 4." << endl; continue; }
            int n;
            cout << "Cuantos patrones quiere ver (n): ";
            cin >> n;
            if (n < 1) {
                cout << "n invalido." << endl; continue;
            }
            cout << "Secuencia: " << idActual << endl;
            topNKmers(hk, n);

        } else if (opcion == 8) {
            string id;
            double beta;
            cout << "ID de la secuencia: ";
            cin >> id;
            int pos = HashTableBusqueda(ht, id);
            if (pos == -1) { cout << "No encontrada." << endl; continue; }
            cout << "Valor de beta: ";
            cin >> beta;
            matrizAdyacencia(ht.bins[pos].dato.sequence, beta);
        }
    }
    return 0;
}
