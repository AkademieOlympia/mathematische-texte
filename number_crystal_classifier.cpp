#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <map>
#include <algorithm>

using namespace std;

struct Factor {
    long long p;
    int exp;
};

vector<Factor> factorize(long long n) {
    vector<Factor> f;
    if (n <= 1) return f;
    int c = 0;
    while (n % 2 == 0) { c++; n /= 2; }
    if (c) f.push_back({2, c});
    for (long long p = 3; p * p <= n; p += 2) {
        c = 0;
        while (n % p == 0) { c++; n /= p; }
        if (c) f.push_back({p, c});
    }
    if (n > 1) f.push_back({n, 1});
    return f;
}

string factor_string(const vector<Factor>& f) {
    if (f.empty()) return "1";
    stringstream ss;
    for (size_t i = 0; i < f.size(); ++i) {
        if (i) ss << " * ";
        ss << f[i].p;
        if (f[i].exp > 1) ss << "^" << f[i].exp;
    }
    return ss.str();
}

char abce_class(long long p) {
    if (p == 2) return '2';
    if (p == 3) return '3';
    long long r = p % 12;
    if (r == 5) return 'A';
    if (r == 7) return 'B';
    if (r == 11) return 'C';
    if (r == 1) return 'E';
    return '?';
}

string csv_escape(const string& s) {
    bool need = false;
    for (char c : s) if (c == ',' || c == '"' || c == '\n') need = true;
    if (!need) return s;
    string out = "\"";
    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out += c;
    }
    out += "\"";
    return out;
}

struct Classification {
    long long n;
    string factors;
    int shell_sigma = 0;
    int alpha2 = 0;
    int alpha3 = 0;
    int vA = 0;
    int vB = 0;
    int vC = 0;
    int vE = 0;
    int v13 = 0;
    long long smooth_core = 1;
    long long defect = 1;
    long long largest_prime_factor = 1;
    string defect_factorization;
    string defect_type;
    string crystal_type;
};

Classification classify(long long n, long long B = 13) {
    Classification c;
    c.n = n;
    auto f = factorize(n);
    c.factors = factor_string(f);

    vector<Factor> defect_factors;

    for (auto &x : f) {
        c.shell_sigma += x.exp;
        c.largest_prime_factor = max(c.largest_prime_factor, x.p);

        if (x.p == 2) c.alpha2 += x.exp;
        else if (x.p == 3) c.alpha3 += x.exp;
        else {
            char cls = abce_class(x.p);
            if (cls == 'A') c.vA += x.exp;
            else if (cls == 'B') c.vB += x.exp;
            else if (cls == 'C') c.vC += x.exp;
            else if (cls == 'E') c.vE += x.exp;
        }

        if (x.p == 13) c.v13 += x.exp;

        if (x.p <= B) {
            for (int i = 0; i < x.exp; ++i) c.smooth_core *= x.p;
        } else {
            for (int i = 0; i < x.exp; ++i) c.defect *= x.p;
            defect_factors.push_back(x);
        }
    }

    c.defect_factorization = factor_string(defect_factors);

    if (n == 1) {
        c.defect_type = "Einheit";
        c.crystal_type = "Ursprung / leere Schale";
    } else if (c.defect == 1) {
        c.defect_type = "rein glatt";
        c.crystal_type = "perfekter Zahlkristall";
    } else {
        int defect_sigma = 0;
        int distinct_defect_primes = (int)defect_factors.size();
        for (auto &d : defect_factors) defect_sigma += d.exp;

        if (distinct_defect_primes == 1 && defect_sigma == 1) {
            char cls = abce_class(defect_factors[0].p);
            string label;
            if (cls == 'A') label = "A-halbglatt";
            else if (cls == 'B') label = "B-halbglatt";
            else if (cls == 'C') label = "C-halbglatt";
            else if (cls == 'E') label = "E-halbglatt";
            else label = "halbglatt";
            c.defect_type = label;
            c.crystal_type = "glatter Kern + einfacher " + string(1, cls) + "-Defekt";
        } else if (distinct_defect_primes == 1 && defect_sigma > 1) {
            char cls = abce_class(defect_factors[0].p);
            c.defect_type = string(1, cls) + string("-Potenzdefekt");
            c.crystal_type = "glatter Kern + Defektpotenz";
        } else {
            c.defect_type = "mehrfach defekt";
            c.crystal_type = "glatter Kern + Defektcluster";
        }
    }

    return c;
}

void write_csv(const string& path, int max_n, long long B) {
    ofstream out(path);
    out << "n,factorization,shell_sigma,alpha2,alpha3,vA,vB,vC,vE,v13,smooth_bound_B,smooth_core,defect,defect_factorization,largest_prime_factor,defect_type,crystal_type\n";
    for (int n = 1; n <= max_n; ++n) {
        auto c = classify(n, B);
        out << c.n << ","
            << csv_escape(c.factors) << ","
            << c.shell_sigma << ","
            << c.alpha2 << ","
            << c.alpha3 << ","
            << c.vA << ","
            << c.vB << ","
            << c.vC << ","
            << c.vE << ","
            << c.v13 << ","
            << B << ","
            << c.smooth_core << ","
            << c.defect << ","
            << csv_escape(c.defect_factorization) << ","
            << c.largest_prime_factor << ","
            << csv_escape(c.defect_type) << ","
            << csv_escape(c.crystal_type) << "\n";
    }
}

void print_classification(long long n, long long B) {
    auto c = classify(n, B);
    cout << "n = " << c.n << "\n";
    cout << "Faktorisierung: " << c.factors << "\n";
    cout << "Schale sigma(n): " << c.shell_sigma << "\n";
    cout << "Sonderachsen: alpha2=" << c.alpha2 << ", alpha3=" << c.alpha3 << "\n";
    cout << "ABCE-Vektor: A=" << c.vA << ", B=" << c.vB << ", C=" << c.vC << ", E=" << c.vE << "\n";
    cout << "E-Tiefe v13(n): " << c.v13 << "\n";
    cout << "Glaettegrenze B: " << B << "\n";
    cout << "Glatter Kern g_B(n): " << c.smooth_core << "\n";
    cout << "Defekt d_B(n): " << c.defect << " = " << c.defect_factorization << "\n";
    cout << "Groesster Primfaktor P+(n): " << c.largest_prime_factor << "\n";
    cout << "Klassifikation: " << c.defect_type << "\n";
    cout << "Kristalltyp: " << c.crystal_type << "\n";
}

int main(int argc, char* argv[]) {
    long long B = 13;        // kleinste vollstaendige ABCE-Basis: 2,3,5,7,11,13
    int max_n = 1000;
    string csv_path = "number_crystal_1_to_1000.csv";

    if (argc >= 2) {
        long long n = stoll(argv[1]);
        if (argc >= 3) B = stoll(argv[2]);
        print_classification(n, B);
        return 0;
    }

    write_csv(csv_path, max_n, B);
    cout << "CSV erzeugt: " << csv_path << " fuer n=1.." << max_n << " mit B=" << B << "\n";
    cout << "Einzelpruefung: ./number_crystal_classifier <n> [B]\n";
    return 0;
}
