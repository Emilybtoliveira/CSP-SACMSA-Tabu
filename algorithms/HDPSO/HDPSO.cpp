using namespace std;
#include <bits/stdc++.h>
#include <chrono>

#define MAX_DIMENSIONS 10000
#define MAX_PARTICLES 10000


class Swap{
public:
    int column; 
    int new_char;

    Swap() {}
};

class Particle {
public:
    int fitness; 
    int best;
    int position[MAX_DIMENSIONS]; 
    int best_position[MAX_DIMENSIONS]; 
    vector<Swap> velocity; 

    void printBest(int size){
        cout << "best particle position: " << endl;

        for(int i = 0; i < size; i++){
            printf("%d ", best_position[i]);
        }
        cout << endl;
    }

    void printPosition(int size){
        cout << "particle position: " << endl;

        for(int i = 0; i < size; i++){
            printf("%d ", position[i]);
        }
        cout << endl;
    }

    void printVelocity(){
        printf("velocity size %d.\n", (int)velocity.size());
        for (int i = 0; i < (int)velocity.size(); i++)
        {
            printf("column: %d new_char: %d\n", velocity[i].column, velocity[i].new_char);
        }        
    }

    Particle() {}
};


// Instance data structures
int n, m, t, min_alpha, max_alpha; // numero de cadeias, tamanho de cada cadeia, tamanho do alfabeto
vector<char> dataset_alphabet;
vector<string> strings_dataset;
vector<vector<int>> integer_dataset;
vector<vector<int>> char_possibilities_per_column;
map<char, int> alphabetMap;
double alphabet_frequency[MAX_DIMENSIONS][MAX_DIMENSIONS];
int** hamming_matrix;


// Auxiliary data structures
int g_best, 
swarm_size = 60, 
loops_with_no_improval = 0, 
loops = 0;
int g_best_position[MAX_DIMENSIONS];
array<Particle, MAX_PARTICLES> particleSwarm;

random_device rd;
mt19937 gen(rd());
unsigned int seed;

// Hyperparameters
double w = 0.8, 
c1 = 8, 
c2 = 4, 
r1, r2, cognitive_weight, social_weight;

//==================== DEBUG FUNCTIONS ===============//

void printBest(){
    cout << "best particle: " << endl;

    for(int i = 0; i < m; i++){
        printf("%d ", g_best_position[i]);
    }
    cout << endl;
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

    for(int i = 0; i < t; i++) {
        for(int j = 0; j < m; j++){
            alphabet_frequency[i][j] = 0;
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
            alphabet_frequency[element][j] += 1;

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

//==================== AUXILIARY FUNCTIONS ===============//

void UpdateBestHamMatrix()
{
    int evaluation = 0;

    for (int i = 0; i < n; i++)
    {
        hamming_matrix[i][m] = 0;
        
        for (int j = 0; j < m; j++){
            evaluation = (integer_dataset[i][j] != g_best_position[j]);
            hamming_matrix[i][j] = evaluation;
            hamming_matrix[i][m] += evaluation;
        }

    }
} 

void InitializeHamMatrix()
{
    int** matrix = new int*[n];
    
    for (int i = 0; i < n; i++)
    {
        matrix[i] = new int[m+1];

        for (int j = 0; j < m+1; j++)
        {
            matrix[i][j] = 0;
        }
    }
    
    hamming_matrix = matrix;
}

//==================== HDPSO FUNCTIONS ===============//
void CopyParticle(int* originalVector, int* cloneVector){
    for (int j = 0; j < m; j++){
        cloneVector[j] = originalVector[j];
    }
}

void SolveSwapConflict(vector<Swap> &velocity, Swap &new_swap){
    for (int i = 0; i < (int)velocity.size(); i++)
    {
        if(velocity[i].column == new_swap.column){
            double prob = (double) rand()/RAND_MAX;
            
            if(prob > 0.5){ //substitui pelo novo
                velocity.erase(velocity.begin() + i);
                velocity.push_back(new_swap);
            }   

            return;             
        }
    }

    velocity.push_back(new_swap);
} 

void InitParticleSwarm(){
    for(int i = 0; i < swarm_size; i++){
        Particle new_particle;
        
        for(int j = 0; j < m; j++){
            vector<double> lotery(t);

            for (int k = 0; k < t; k++) {
                lotery[k] = alphabet_frequency[k][j] / (double) m;
            }

            discrete_distribution<int> column_distribution(lotery.begin(), lotery.end());


            int character = column_distribution(gen);
            new_particle.position[j] = character;
            new_particle.best_position[j] = character;
            new_particle.best = m;
        }

        particleSwarm[i] = new_particle;
    }

    // for(int i = 0; i < swarm_size; i++){
    //     Particle new_particle = particleSwarm[i];
    //     new_particle.printPosition(m);
    // }
}

int EvaluateFitness(int* particle){
    int distance = 0, max_hamming = 0;

    for (int i = 0; i < n; i++)
    {
        distance = 0;

        for (int j = 0; j < m; j++)
            distance += (integer_dataset[i][j] != particle[j]);        

        if (distance > max_hamming)
           max_hamming = distance;
    }
    

    return max_hamming;
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

void BestImprovementLS(){
    int best = m, opt, current_char, best_char = 0, best_column = 0;
    bool improval = true;

    while(improval){
        best = m;

        for(int i = 0; i < m; i++){            
            for(int k = 0; k < (int)char_possibilities_per_column[i].size(); k++){
                if(g_best_position[i] != char_possibilities_per_column[i][k]){
                    current_char = g_best_position[i];
                    g_best_position[i] = char_possibilities_per_column[i][k];

                    opt = CalculateHammingAtNeighbourhood(i, current_char, g_best_position[i]);  

                    if(opt < best){
                        //printf("opt %d  best %d\n", opt, best);
                        best_column = i;
                        best_char = char_possibilities_per_column[i][k];
                        best = opt;
                    }

                    g_best_position[i] = current_char;                
                }
            }
        }

        if(best < g_best){
            //printf("best %d  g_best %d\n", best, g_best);
            g_best_position[best_column] = best_char;
            g_best = best;
            return;
        } else {
            improval = false;
        }
    }
}

vector<Swap> GetSwapArray(int* particle_position, int* best_particle_position){
    vector<Swap> swapArray;

    for(int i = 0; i < m; i++){
        if(particle_position[i] != best_particle_position[i]){
            Swap new_swap;
            new_swap.column = i;
            new_swap.new_char = best_particle_position[i];

            swapArray.push_back(new_swap);
        }
    }

    return swapArray;
}

void CalculateVelocity(Particle &particle){
    vector<Swap> new_velocity;
    vector<Swap> cognitive_swap = GetSwapArray(particle.position, particle.best_position);
    vector<Swap> social_swap = GetSwapArray(particle.position, g_best_position);

    if(particle.velocity.size()) { /* componente inercial */
        double inertial_weight = round(w * particle.velocity.size());
        vector<int> inertial_indexes;
        
        for(int i = 0; i < (int)particle.velocity.size(); i++){
            inertial_indexes.push_back(i);
        }

        shuffle(inertial_indexes.begin(), inertial_indexes.end(), default_random_engine(seed));

        for(int i = 0; i < inertial_weight; i++){
            new_velocity.push_back(particle.velocity[inertial_indexes[i]]);
        }
    }

    if(cognitive_swap.size() && cognitive_weight){ /* componente cognitivo */
        vector<int> cognitive_indexes;
        
        for(int i = 0; i < (int)cognitive_swap.size(); i++){
            cognitive_indexes.push_back(i);
        }

        shuffle(cognitive_indexes.begin(), cognitive_indexes.end(), default_random_engine(seed));

        for(int i = 0; i < min((int)cognitive_weight, (int)cognitive_swap.size()); i++){
            SolveSwapConflict(new_velocity, cognitive_swap[cognitive_indexes[i]]);
            //new_velocity.push_back(cognitive_swap[cognitive_indexes[i]]);
        }
    }

    if(social_swap.size() && social_weight){ /* componente social */
        vector<int> social_indexes;

        for(int i = 0; i < (int)social_swap.size(); i++){
            social_indexes.push_back(i);
        }

        shuffle(social_indexes.begin(), social_indexes.end(), default_random_engine(seed));
    
        for(int i = 0; i < min((int)social_weight, (int)social_swap.size()); i++){
            SolveSwapConflict(new_velocity, social_swap[social_indexes[i]]);
            //new_velocity.push_back(social_swap[social_indexes[i]]);
        }
    }      

    particle.velocity = new_velocity;
}

void MoveParticle(Particle &particle){
    Swap swap;
    
    CalculateVelocity(particle);
    
    for(int i = 0; i < (int) particle.velocity.size(); i++){
        swap = particle.velocity[i];

        particle.position[swap.column] = swap.new_char;
    }
}

//==================== CORE ===============//

void HDPSO()
{
    auto loop_start = chrono::high_resolution_clock::now();
    auto loop_end = chrono::high_resolution_clock::now();
    auto loop_cur = chrono::duration_cast<chrono::milliseconds>(loop_end - loop_start);  

    bool improve = false;
    g_best = m;

    InitializeHamMatrix();
    InitParticleSwarm();

    while(loop_cur.count() < 30000)
    {
        improve = false;
        seed = chrono::system_clock::now().time_since_epoch().count();
        srand(seed);
        
        r1 = (double) rand() / RAND_MAX;
        r2 = (double) rand() / RAND_MAX;        

        cognitive_weight = round(c1 *r1);
        social_weight = round(c2 * r2);
        
        for(int i = 0; i < swarm_size; i++){
            particleSwarm[i].fitness = EvaluateFitness(particleSwarm[i].position);
            
            if(particleSwarm[i].fitness < particleSwarm[i].best){
                //atualiza o best_particle
                CopyParticle(particleSwarm[i].position, particleSwarm[i].best_position);
                particleSwarm[i].best = particleSwarm[i].fitness;                

                if(particleSwarm[i].fitness < g_best){
                    //atualiza o g_best_position

                    CopyParticle(particleSwarm[i].position, g_best_position);
                    g_best = particleSwarm[i].fitness;

                    UpdateBestHamMatrix();   
                    //PrintHamMatrix();                 
                    BestImprovementLS();
                    
                    improve = true;
                    loops_with_no_improval = 0;
                }
            }
            
            MoveParticle(particleSwarm[i]);
        }

        if(!improve) loops_with_no_improval++;        
        loops++;

        loop_end = chrono::high_resolution_clock::now();
        loop_cur = chrono::duration_cast<chrono::milliseconds>(loop_end - loop_start);      
    }


    cout << "opt: " << g_best << endl;
    cout << "time: " << loop_cur.count() << endl;
    cout << "loops: " << loops << endl;
}

int main()
{
    ReadInstance();
    GenerateAlphabetMapping();
    InstanceTransformFunc();
    //PrintInstance();
    
    HDPSO();   


    return 0;
}