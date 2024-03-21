#include <bits/stdc++.h>
#include <chrono>
using namespace std;

int n, m, t, min_alpha, max_alpha; // numero de cadeias, tamanho de cada cadeia, tamanho do alfabeto
vector<char> dataset_alphabet;
vector<string> strings_dataset;
vector<vector<int>> integer_dataset;
map<char, int> alphabetMap;
unsigned seed;

vector<int> farthest_string;
vector<int> best_t;
vector<int> most_frequent;

int time_limit;
int delta, best_t_energy, current_t_energy, loops = 0;    


//==================== DEBUG FUNCTIONS ===============//

void printH(int **H, int j_limit){
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= j_limit; j++)
        {
            cout << H[i][j] << " ";
        }
        cout << endl;        
    }    
}

void printT(vector<int> t){
    for(int i = 0; i < m; i++) cout << t[i] << " ";
    cout << endl;
}

//==================== INSTANCE FUNCTIONS ===============//

void PrintInstance()
{
    printf("%distance %distance %distance\n", n, m, t);

    for (char i : dataset_alphabet)
    {
        cout << i << '\n';
    }

    for (const auto &item : alphabetMap)
    {
        cout << "[" << item.first << ", " << item.second << "]\n";
    }

    for (string i : strings_dataset)
    {
        cout << i << '\n';
    }

    for (vector<int> integer : integer_dataset)
    {
        for(int j = 0; j < m; j++){
            cout << integer[j] << ' ';
        }
        cout << endl;
    }    
}

void GenerateAlphabetMapping()
{
    int i;
    for (i = 0; i < t; i++)
    {  
        char alpha_char = dataset_alphabet[i];
        alphabetMap[alpha_char] = i;
    }

    min_alpha = 0;
    max_alpha = i - 1;
}

void InstanceTransformFunc()
{
    vector<int> int_string;
    int element;


    for (string cur_string : strings_dataset)
    {
        int_string.clear();

        for (int j = 0; j < (int)cur_string.size(); j++)
        {
            element = alphabetMap[cur_string[j]];
            int_string.push_back(element);            
        }
        integer_dataset.push_back(int_string);
    }   
}

void ReadInstance()
{
    //cin >> t >> n >> m ; // tamanho do alfabeto, numero de cadeias, tamanho das cadeias
    cin  >> n >> m >> t; //  numero de cadeias, tamanho das cadeias, tamanho do alfabeto

    char cur_char;

    for (int i = 0; i < t; i++)
    {   
        cin >> cur_char;
        dataset_alphabet.push_back(cur_char);
    }

    for (int i = 0; i < n; i++)
    {
        string cur_string;
        cin >> cur_string;
        strings_dataset.push_back(cur_string);
    }
}

//==================== SAGW FUNCTIONS ===============//

int JthMostFrequent(int j){
    int count[t];
    int prob, most_freq = 0;

    for (int i = 0; i < t; i++) count[i] = 0;
    for (int i = 0; i < n; i++) count[integer_dataset[i][j]] += 1;
    
    for (int i = 0; i < t; i++)
    {
        if(count[i] > count[most_freq])
        {
            most_freq = i;
        }
        else if(count[i] == count[most_freq])
        {
            prob = rand() % 2; 
            if(prob) most_freq = i; 
        } 
    }       

    return most_freq;
}

int ChooseJthChar(int j, int **H){       
    if(j == 0) //return most frequent char
    {   
        return JthMostFrequent(j);
    } else //greedy walk
    {
        int max_ham = 0;
        vector<int> max_ham_positions;

        for (int i = 0; i < n; i++)
        {
            if (H[i][j-1] > max_ham) max_ham = H[i][j-1];
        }

        for (int i = 0; i < n; i++)
        {
            if (H[i][j-1] == max_ham) max_ham_positions.push_back(i);
        }

        int pos = rand() % max_ham_positions.size();
        int row = max_ham_positions[pos];
        int character = integer_dataset[row][j]; 
        
        //printf("%d, %d, %d, %d\n", max_ham, pos, row, character);

        return character;
    }
}

vector<int> GreedyWalk(){
    vector<int> t;
    
    int **H = new int *[n];
    for(int i = 0; i < n; i++)
        H[i] = new int[m];
        
    for (int j = 0; j < m; j++)
    {
        t.push_back(ChooseJthChar(j, H));

        for (int i = 0; i < n; i++){
            if(j > 0) H[i][j] = (integer_dataset[i][j] != t[j]) + H[i][j-1];
            else  H[i][j] = (integer_dataset[i][j] != t[j]);
        }
        //printH(H, j);
    }

    return t; 
}

vector<int> Mutate(vector<int> t, float mutation_rate){
    float prob;
    int most_freq;
    int flag = 0;

    for (int j = 0; j < m; j++) // tenta substituir pelo char mais frequente
    {
        prob = (float) rand() / RAND_MAX;

        if (prob < mutation_rate) 
        {
            most_freq = JthMostFrequent(j);
            if(most_freq != t[j]) flag = 1; 
            t[j] = most_freq;
        } 
    }

    if(!flag) // se nao der certo, substitui pela farthest string
    {
        //cout << "No mutation, using farthest string.\n";

        for (int j = 0; j < m; j++) 
        {
            prob = (float) rand() / RAND_MAX;

            if (prob < mutation_rate) 
            {   
                t[j] = farthest_string[j];
            } 
        }
    } 

    return t;

}

int EvaluateEnergy(vector<int> t){
    int distance = 0, max_hamming = 0;

    for (vector<int> S : integer_dataset)
    {
        distance = 0;

        for (int i = 0; i < m; i++)
        {
            if (S[i] != t[i])
            {
                distance++;
            }
        }

        if (distance > max_hamming)
        {
            max_hamming = distance;
            farthest_string = S;
        }
    }

    return max_hamming;
}

//==================== CORE ===============//

void compare(){
    for(int i = 0; i < m; i++){
        most_frequent.push_back(JthMostFrequent(i));
    }
    int count = 0;

    for(int i = 0; i < m; i++){
        if(most_frequent[i] == best_t[i]) count++;
    }

    printf("using most frequent in %d/%d\n", count, m);
}

void SAGW(int T_max, int L, float gama, float mutation_rate){
    float prob;
    vector<int> current_t;
    float T = T_max;
    int loops_with_no_improval = 0;   

    auto loop_start = chrono::high_resolution_clock::now();
    auto loop_end = chrono::high_resolution_clock::now();
    auto loop_cur = chrono::duration_cast<chrono::milliseconds>(loop_end - loop_start);  


    best_t = GreedyWalk();
    best_t_energy = EvaluateEnergy(best_t);

    //cout << best_t_energy << endl;

    while(loop_cur.count() < 30000)
    {    
        for (int i = 0; i < L; i++) 
        {
            seed = chrono::system_clock::now().time_since_epoch().count(); 
            srand(seed); 

            current_t = Mutate(best_t, mutation_rate);
            current_t_energy = EvaluateEnergy(current_t);

            
            delta =  current_t_energy - best_t_energy;
            prob = (float) rand() / RAND_MAX;


            if((delta <= 0) or (exp(-delta/T) > prob))
            {
                best_t = current_t;
                best_t_energy = current_t_energy;
                loops_with_no_improval = 0;
            } else {
                loops_with_no_improval += 1;
            }
        }
        
        loops++;
        T = T * gama;

        loop_end = chrono::high_resolution_clock::now();
        loop_cur = chrono::duration_cast<chrono::milliseconds>(loop_end - loop_start);  
    } 

        
    cout << "opt: " << best_t_energy << endl;
    cout << "time: " << loop_cur.count() << endl;
    cout << "loops: " << loops*500 << endl; 
}

int main(int argc, char *argv[])
{
    for (int i = 0; i < argc; i++) 
    {        
        if(!strcmp(argv[i], "-t")) {          
            sscanf(argv[i+1], "%d", &time_limit);
            i++;
        }     
    }

    ReadInstance();
    GenerateAlphabetMapping();
    InstanceTransformFunc();
    //PrintInstance();

    SAGW(5, 500, 0.8, 0.4);
    

    return 0;
}
