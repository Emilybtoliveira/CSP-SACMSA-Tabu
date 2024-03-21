#include <ilcplex/ilocplex.h>
#include <bits/stdc++.h>
#include <chrono>

using namespace std;
#define MAX_SETS 100000 //NOTE: Limita o tamanho das strings em MAX_SETS/(max_age+1) = 50000
#define MAX_DIMENSIONS 10000
#define MAX_SET_SIZE 1000
#define MAX_ALPHABET 26
#define INF -574873548

//========================== data structures ===========================

class Set {
public:
    int column;
    int closest_string;
    int selected;
    array<int, MAX_DIMENSIONS> hamming_dist; 

    Set() : selected(0) {}
};

class HashMap {
public:
    void insert(int key, const Set& set) {
        table_[key] = set;
    }

    int size() {
        return table_.size();
    }

    Set* find(int key) {
        auto iter = table_.find(key);
        if (iter != table_.end()) {
            return &(iter->second);
        }
        return nullptr;
    }

    void remove(int key) {
        table_.erase(key);
    }

public:
    unordered_map<int, Set> table_;
};


/* INSTANCE PARAMETERS */
int n, m, t; // numero de cadeias, tamanho de cada cadeia, tamanho do alfabeto
vector<char> dataset_alphabet;
vector<string> strings_dataset;
vector<vector<int>> integer_dataset;
map<char, int> alphabetMap;
vector<vector<int>> char_possibilities_per_column;
vector<int> frequency_per_column;

/*AUXILIARY DATASTRUCTURES */
HashMap columns_sets;
int current_index = 0,
loops_with_no_improval = 0,
loops = 0,
n_new_sets = 0,
n_a = 1;

double bsf, opt, prob;
int character, key, column, known_opt = 0;

int** hamming_matrix;
int visited[MAX_ALPHABET][MAX_SETS]; 
int tabu[MAX_ALPHABET][MAX_SETS]; 
int mapping[MAX_SETS];
int s_bsf[MAX_DIMENSIONS]; 

/* HYPERPARAMETERS */
int tabu_threshold = 0;
float t_ILP = 3,    /* Tempo limite pro exato */
t_prop = 0.5596, /* Proporção aceitável de tempo pra resolver */
t_solve = 0, /* Tempo que o exato levou para resolver (não é hiperparametro) */
alpha_bsf = 0, /* Taxa de geração de componentes (não é hiperparametro) */
alpha_UB = 0.9872, /* Valor máximo do alpha_bsf */
alpha_LB = 0.2950, /* Valor minimo do alpha_bsf */
alpha_red = 0.9385; /* Fator de redução do alpha_bsf */

//==================== DEBUG FUNCTIONS ===============//

void print_hash() {
    cout << "Columns_sets size: " << columns_sets.size() << endl;

    for (const auto& pair : columns_sets.table_) {
        cout << "Key: " << pair.first << endl;
        const Set& set = pair.second;
        
        cout << "column: " << set.column << "\n";
        
        cout << "hamming distances: [ ";
        for (int i = 0; i < n; i++) {
            cout << set.hamming_dist[i] << " ";
        }
        cout << "]" << endl;
        
        cout << "closest string: " << set.closest_string << "\n";

        cout << endl;
    }
}

void activeChars(){
    cout << "[";

    for(int j = 0; j < m; j ++){
        cout << visited[t][j] << ", ";
    }

    cout << "]\n";
}

void PrintHamMatrix()
{
   for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++){
            printf("%d ", hamming_matrix[i][j]);
        }

        printf("sum: %d\n", hamming_matrix[i][m]);
    } 
}

//==================== INSTANCE FUNCTIONS ===============//

void generateAlphabetMapping()
{
    int i;
    for (i = 0; i < t; i++)
    {
        char alpha_char = dataset_alphabet[i];
        alphabetMap[alpha_char] = i;
    }

    /* for (const auto &item : alphabetMap)
    {
        cout << "[" << item.first << ", " << item.second << "]\n";
    } */
}

int instanceTransformFunc()
{
    int frequency[t][m];
    vector<int> int_string;
    int element, max_freq, max_freq_index, smallest_char_per_col = t;

    for(int i = 0; i < t; i++) {
        for(int j = 0; j < m; j++){
            frequency[i][j] = 0;
        }
    }

    for (int j = 0; j < m; j++){
        char_possibilities_per_column.push_back(int_string);
    }

    for (string cur_string : strings_dataset)
    {
        int_string.clear();

        for (int j = 0; j < (int)cur_string.size(); j++)
        {
            element = alphabetMap[cur_string[j]];
            int_string.push_back(element);
            frequency[element][j] += 1;

            if(char_possibilities_per_column[j].size() > 0){                
                if (find(char_possibilities_per_column[j].begin(), char_possibilities_per_column[j].end(), element) == char_possibilities_per_column[j].end()){
                    char_possibilities_per_column[j].push_back(element);  
                } 
            }
            else{
                char_possibilities_per_column[j].push_back(element);  
            }
        }

        integer_dataset.push_back(int_string);
    }
    

    for(int j = 0; j < m; j ++) { /* Counts most frequent per column*/
        max_freq = 0;
        max_freq_index = 0;

        prob = rand() % 2;

        for(int i = 0; i < t; i++){
            if (frequency[i][j] > max_freq || (frequency[i][j] == max_freq && prob) ){
                max_freq = frequency[i][j];
                max_freq_index = i;
            } 
        }

        frequency_per_column.push_back(max_freq_index);
    }

    for (int i = 0; i < m; i++){      /* Calculates de column that has least options */   
        if(char_possibilities_per_column[i].size() < smallest_char_per_col){
            smallest_char_per_col = char_possibilities_per_column[i].size(); 
        }
    }

    return smallest_char_per_col;
}

void readInstance()
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

//==================== AUXILIARY FUNCTIONS ===============//

void initializeDataStructures(){
    hamming_matrix = new int*[n];
    
    for (int i = 0; i < n; i++)
    {
        hamming_matrix[i] = new int[m+1];

        for (int j = 0; j < m+1; j++)
        {
            hamming_matrix[i][j] = 0;
        }
    }
    
    
    for(int i = 0; i < MAX_ALPHABET;i++){
        for(int j = 0; j < MAX_SETS; j++){
            visited[i][j] = 0;
            tabu[i][j] = INF;
        }
    }

    for(int i = 0; i < m; i++){
        tabu[t][i] = 0;
    }

    srand(chrono::system_clock::now().time_since_epoch().count());
}

void updateTabuList(){
    for(int i = 0; i < t; i++){
        for(int j = 0; j < m; j++){
            if((tabu[i][j]!= INF) && (loops - tabu[i][j] > tabu_threshold)){
                tabu[i][j] = INF;
                tabu[t][j] -= 1;
            }
        }
    }
}

int getRandomChar(int column){
    int sorted_pos = rand() % (char_possibilities_per_column[column].size());  //escolhe aleatoriamente um caractere para a posição
    character = char_possibilities_per_column[column][sorted_pos];

    while ( (visited[character][column] == 1) or (loops - tabu[character][column] <= tabu_threshold)){ 
        /* aqui o loop varre na direção antihorario as linhas de caracteres até encontrar uma posição em que o caractere n esteja sendo utilizado
            e que nao tenha sido utilizado ha pelo menos tabu_threshold loops */      
        sorted_pos -= 1;
        if (sorted_pos == -1) sorted_pos = char_possibilities_per_column[column].size() -1; 
        character = char_possibilities_per_column[column][sorted_pos];
    }

    return character;
}

void addHammingToSet(Set* set, int character, int column, int flag){
    for(int i = 0; i < n; i++){
        set->hamming_dist[i] = (character != integer_dataset[i][column]);

        if(flag){
            hamming_matrix[i][column] = set->hamming_dist[i];
            hamming_matrix[i][m] += set->hamming_dist[i];
        }
    }
}

void UpdateHamMatrix(array<int, MAX_DIMENSIONS> ham_distances, int column)
{
    for (int i = 0; i < n; i++)
    {
        if(!hamming_matrix[i][column] && ham_distances[i]) 
            hamming_matrix[i][m] += 1;
        if(hamming_matrix[i][column] && !ham_distances[i]) 
            hamming_matrix[i][m] -= 1;    
        
        hamming_matrix[i][column] = ham_distances[i];
    }
} 

int CalculateHammingAtNeighbourhood(int column, int previous_char, int new_char){
    vector<int> hammings;
    for(int i = 0; i < n; i++){
        hammings.push_back(hamming_matrix[i][m]);
    }
    
    for(int i = 0; i < n; i++){
        if(integer_dataset[i][column] != previous_char && integer_dataset[i][column] == new_char)
            hammings[i] -= 1;
        else if(integer_dataset[i][column] == previous_char && integer_dataset[i][column] == new_char)
            hammings[i] += 1;
    }

    return *max_element(hammings.begin(), hammings.end());;
}

//===================SAGW FUNCTIONS =======================//

int evalHamDistance(int* s){
    int distance = 0, max_hamming = 0;

    for (int i = 0; i < n; i++)
    {
        distance = 0;

        for (int j = 0; j < m; j++)
        {
            if (integer_dataset[i][j] != s[j]) distance++;
        }

        if (distance > max_hamming)
        {
            max_hamming = distance;
        }
    }

    return max_hamming;
}

int chooseJthChar(int j){
    if(!visited[frequency_per_column[j]][j]) 
        return frequency_per_column[j]; 

    return getRandomChar(j);
}

int computeJthChar(int j){       
    if(j == 0) //return most frequent char
    {   
        return frequency_per_column[j];
    } else //greedy walk
    {
        int max_ham = 0;
        vector<int> max_ham_positions;

        for (int i = 0; i < n; i++)
        {
            if (hamming_matrix[i][j-1] > max_ham) max_ham = hamming_matrix[i][j-1];
        }

        for (int i = 0; i < n; i++)
        {
            if (hamming_matrix[i][j-1] == max_ham) max_ham_positions.push_back(i);
        }

        character = rand() % max_ham_positions.size();
        
        return integer_dataset[max_ham_positions[character]][j];
    }
}

/* SAGW initial solution */
void GWInitialSolution(){           
    for (int j = 0; j < m; j++)
    {
        Set new_data;
        columns_sets.insert(j, new_data);
    }

    for (int j = 0; j < m; j++)
    {
        Set* new_data = columns_sets.find(j);
        character = computeJthChar(j);
        s_bsf[j] = character;

        new_data->column = j;
        new_data->closest_string = character;        
        new_data->selected = 1;

        visited[character][j] = 1;
        visited[t][j] += 1;
        
        for (int i = 0; i < n; i++){
            if(j > 0) hamming_matrix[i][j] = (integer_dataset[i][j] != character) + hamming_matrix[i][j-1];
            else  hamming_matrix[i][j] = (integer_dataset[i][j] != character);
        }
    }

    current_index += m;
}

/* Most frequent initial solution */
void MFInitialSolution(){    
    for(int j = 0; j < m; j++){
        Set new_data;
        columns_sets.insert(j, new_data);
    }

    for(int j = 0; j < m; j++){
        Set* new_data = columns_sets.find(j);
        character = frequency_per_column[j]; //usando a heuristica gulosa do mais frequente
        s_bsf[j] = character;

        new_data->column = j;
        new_data->closest_string = character;
        new_data->selected = 1;

        addHammingToSet(new_data, character, j, 1);

        visited[character][j] = 1;
        visited[t][j] += 1;    
    }

    current_index += m;
}


//==================== SELF ADAPT CMSA FUNCTIONS ===============//
/* ilomipex4: Spend at most timeLimit sec. on optimization OR quit as soon as an acceptable (better) solution is found 
ILOMIPINFOCALLBACK: https://www-eio.upc.es/lceio/manuals/cplex-11/html/refcppcplex/html/macros/ILOMIPINFOCALLBACK0.html */
ILOMIPINFOCALLBACK3(timeLimitCallback,
                    IloCplex, cplex,
                    IloBool,  aborted,
                    IloNum,   timeStart)
{
       
    if ( !aborted  &&  hasIncumbent() ) {
        IloNum timeUsed = cplex.getCplexTime() - timeStart;
        IloNum cost = getIncumbentObjValue();
        
        //printf("Found solution value: %f\n", cost);

        if ( timeUsed >= t_ILP || cost < bsf ) {
            //printf("Cost is better than bsf and time limit run out. Quiting after %fs...\n", timeUsed);
            aborted = IloTrue;
            abort();
        }
    }
}

void bestImprovementLS(){
    int best_n, previous_char, current_n, best_column, best_char;

    while(true){
        best_n = m;
        
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < (int) char_possibilities_per_column[i].size(); j++)
            {
                if(s_bsf[i] != char_possibilities_per_column[i][j]){                    
                    current_n = CalculateHammingAtNeighbourhood(i, s_bsf[i], char_possibilities_per_column[i][j]);
                    
                    if(current_n < best_n){
                        best_column = i;
                        best_char = char_possibilities_per_column[i][j];
                        best_n = current_n;
                    } 
                }
            }        
        }

        if(best_n < bsf){
            s_bsf[best_column] = best_char;
            bsf = best_n;
            cout << bsf << endl;
        } else {
            break;
        }    
    }
}

void updateSelfAdaptParams()
{ 
    if(t_solve < (t_prop*t_ILP) && alpha_bsf > alpha_LB)
        alpha_bsf -= alpha_red;
    if(opt < bsf)
        n_a = 1;
    else {
        if(opt > bsf){
            if(n_a == 1) 
                alpha_bsf = min(alpha_bsf+(alpha_red/10), alpha_UB);
            else
                n_a = 1;
        } else {
            n_a++;
        }
    }
}

void constructAndMerge() 
{
    n_new_sets = 0;

    for (int j = 0; j < m; j++) {
        prob = (float) rand() / RAND_MAX;

        if (prob > alpha_bsf and (visited[t][j] + tabu[t][j]) < (int) char_possibilities_per_column[j].size()) {
            Set new_data;
            new_data.column = j;
            columns_sets.insert(current_index + n_new_sets, new_data);
            n_new_sets++;
        }
    }

    for (int i = current_index; i < current_index + n_new_sets; i++){
        Set* found_data = columns_sets.find(i);
        character = chooseJthChar(found_data->column);

        found_data->closest_string = character;

        addHammingToSet(found_data, character, found_data->column, 0);

        visited[character][found_data->column] = 1;
        visited[t][found_data->column] += 1;
    }

    current_index += n_new_sets;
}

double solve(int n_sets, double bsf) 
{
    IloEnv env;
    IloModel Model(env, "Problema da Seleção de Conjuntos");
    double cost;
    Set* found_data;

    try
    {
        // Variável de decisão
        IloIntVarArray x(env, n_sets, 0, 1);

        // Variável objetivo
        IloNumVar z(env, 0, IloInfinity, ILOINT);

        IloExprArray column;
        column = IloExprArray(env, m);

        for (int j = 0; j < m; j++) column[j] = IloExpr(env);

        int index = 0;
        for (const auto& pair : columns_sets.table_) {
            const Set& set = pair.second;
           
            column[set.column] += x[index];

            mapping[index] = pair.first;
            index += 1;
        }

        for (int j = 0; j < m; ++j) Model.add(column[j] == 1);

        for (int i = 0; i < n; i++)
        {
            IloExpr exp2(env);

            for (int s = 0; s < index; s++)
            {
                key = mapping[s];
                
                found_data = columns_sets.find(key);

                exp2 += found_data->hamming_dist[i] * x[s];
            }

            Model.add(exp2 <= z);
        }

        // restrição p podar mais rapido
        //Model.add(z < bsf);

        // Função objetivo.
        Model.add(IloMinimize(env, z));

        // Solving
        IloCplex cplex(Model);
        //cout << cplex.getModel() << endl;
        cplex.setOut(env.getNullStream());
        //cplex.setParam(IloCplex::Param::TimeLimit, t_ILP); 
        cplex.use(timeLimitCallback(env, cplex, IloFalse, cplex.getCplexTime()));
        //cplex.setParam(IloCplex::EpGap,  0.05); //set the minimum required gap
        //cout << cplex.getModel() << endl;
        
        IloNum time_start = cplex.getCplexTime();

        if (cplex.solve()){
            cost = cplex.getObjValue();
            t_solve = cplex.getCplexTime() - time_start; 
            
            //cout << "solved in time " << t_solve << endl;
            // cout << cplex.getValue(z) << endl ;
            //cout << cplex.getObjValue() << endl << endl;
            //cout << cplex.getCplexStatus() << endl;

            // Obtendo a solução
            IloNumArray sol(env, n_sets);
            cplex.getValues(sol, x);

            for (int i = 0; i < n_sets; i++)
            {
                key = mapping[i];
                found_data = columns_sets.find(key);
                
                if (sol[i] > 0.5){
                    found_data->selected = 1;
                } else{
                    found_data->selected = 0;
                } 
                
            }

            cplex.end();
            Model.end();
            env.end();
       
            return cost;
        }
    }
    catch (const IloException &e)
    {
        cerr << "Exception caught: " << e << endl;
    }
    catch (...)
    {
        cerr << "Unknown exception caught!" << endl;
    }

    return bsf; // se nao encontrou, retorna a melhor que já tem
}

void adapt()
{
    Set* found_data;
    int size = columns_sets.size();

    for (int i = 0; i < size; i++)
    {
        key = mapping[i];
        found_data = columns_sets.find(key);
        character = found_data->closest_string;
        column = found_data->column;

        if (!found_data->selected)
        {  
            visited[character][column] = 0;
            visited[t][column] -= 1;

            tabu[character][column] = loops;  /* se ultrapassou o age, sai do conjunto de soluções e entra na lista tabu */                    
            tabu[t][column] += 1; //ultima linha de cada coluna diz quantos caracteres estão na lista tabu
            //cout << tabu[t][column] << endl;

            columns_sets.remove(key);                    
        } else {
            s_bsf[column] = character;

            if(tabu[character][column] != INF){
                tabu[character][column] = INF; /* se passou a fazer parte da solução, nao faz parte mais da lista tabu */
                tabu[t][column] -= 1; //ultima linha de cada coluna diz quantos caracteres estão na lista tabu
            }

            //UpdateHamMatrix(found_data->hamming_dist, column);
            //PrintHamMatrix();
        }
    }
}

//==================== CORE ===============//

void CMSA()
{
    auto loop_start = chrono::high_resolution_clock::now();
    auto loop_end = chrono::high_resolution_clock::now();
    auto loop_cur = chrono::duration_cast<chrono::milliseconds>(loop_end - loop_start);
    
    alpha_bsf = alpha_UB;
    bsf = m;

    MFInitialSolution();    
    //GWInitialSolution(); 
    
    //bsf = evalHamDistance(s_bsf); 
    //printf("Solução inicial: %d\n", bsf);   
    

    //while(bsf > known_opt and loop_cur.count() < 60000)
    // while (loops_with_no_improval < 5)
    while(loop_cur.count() < 30000)
    {       
        // printf("alpha_bsf: %f\n", alpha_bsf);
        // printf("n_a: %d\n", n_a\);

        for (int i = 0; i < n_a; i++)
            constructAndMerge();           
        
        //cout << "Total de componentes ativos: " << columns_sets.size() << endl;

        opt = solve(columns_sets.size(), bsf);
        //printf("%f, %f\n", bsf, opt);        
        
        updateSelfAdaptParams();
        adapt(); 
        updateTabuList();

        if (opt < bsf)
        {         
            loops_with_no_improval = 0;
            bsf = opt;
            //bestImprovementLS();       
        }
        else
        {
            loops_with_no_improval += 1;
        }    

        srand(chrono::system_clock::now().time_since_epoch().count());
        loops += 1;

        loop_end = chrono::high_resolution_clock::now();
        loop_cur = chrono::duration_cast<chrono::milliseconds>(loop_end - loop_start);              
    }

    cout << "loops: " << loops << endl;
    cout << "time: " << loop_cur.count() << endl;
    cout << "opt: " << bsf << endl;
}

int main(int argc, char *argv[])
{
    readInstance();
    //printInstance();
    generateAlphabetMapping();
    tabu_threshold = instanceTransformFunc();
    //printf("tabu threshold: %d\n", tabu_threshold);    

    initializeDataStructures();    

    //lendo argumentos da linha de comando
    for (int i = 0; i < argc; i++) 
    {               
        if(!strcmp(argv[i], "-tabu")) {          
            sscanf(argv[i+1], "%d", &tabu_threshold);
            //printf("tabu threshold: %d\n", tabu_threshold);
            i++;
        }
        if(!strcmp(argv[i], "-opt")) {          
            sscanf(argv[i+1], "%d", &known_opt);
            //printf("known opt: %d\n", known_opt);
            i++;
        }
        if(!strcmp(argv[i], "-tILP")) {          
            sscanf(argv[i+1], "%f", &t_ILP);
            //printf("t_ILP: %f\n", t_ILP);
            i++;
        }
        if(!strcmp(argv[i], "-tPROP")) {          
            sscanf(argv[i+1], "%f", &t_prop);
            //printf("t_PROP: %f\n", t_prop);
            i++;
        }
        if(!strcmp(argv[i], "-aLB")) {          
            sscanf(argv[i+1], "%f", &alpha_LB);
            //printf("alpha_LB: %f\n", alpha_LB);
            i++;
        }
        if(!strcmp(argv[i], "-aUB")) {          
            sscanf(argv[i+1], "%f", &alpha_UB);
            //printf("alpha_UB: %f\n", alpha_UB);
            i++;
        }
        if(!strcmp(argv[i], "-aRED")) {          
            sscanf(argv[i+1], "%f", &alpha_red);
            //printf("alpha_red: %f\n", alpha_red);
            i++;
        }
    }
    
    CMSA();


    return 0;
}