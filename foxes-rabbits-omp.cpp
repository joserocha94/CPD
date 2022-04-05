#include <vector>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

enum item {FOX, RABBIT, ROCK, NONE};
enum move {UP, RIGHT, DOWN, LEFT};

struct Resident
{
    item type;
    int starving_age;
    int breeding_age;

    Resident() : type(NONE), starving_age(0), breeding_age(0) {}
};

struct Cell 
{
    Resident resident;
    int flag;

    Cell() : resident(), flag(0) {}
};

std::vector<std::vector<Cell>> world;
std::vector<std::vector<Cell>> world_bck;

int N_ROCKS, N_FOXES, N_RABBITS;
int FOX_B_AGE, FOX_S_AGE, RABBIT_B_AGE;
int GENERATIONS;

int FOXES, RABBITS, ROCKS;
int M, N;

float r4_uni(uint32_t *seed)
{
    int seed_input, sseed;
    seed_input = *seed;

    *seed ^= (*seed << 13);
    *seed ^= (*seed >> 17);
    *seed ^= (*seed << 5);
    sseed = *seed;

    return 0.5 + 0.2328306e-09 * (seed_input + sseed);
}

void generate_element(int n, uint32_t *seed, item new_resident, int breeding_age, int starving_age)
{
    int i, j, k;
    for (k = 0; k < n; k++) 
    {
        i = M * r4_uni(seed);
        j = N * r4_uni(seed);
        if (world[i][j].resident.type == NONE)
        {
            world[i][j].resident.type = new_resident;
            world[i][j].resident.breeding_age = breeding_age;
            world[i][j].resident.starving_age = starving_age;
        }
    }
}

void init_world(std::vector<std::vector<Cell>> &cells, long lines, long cols)
{
	cells.clear();
	cells.reserve(lines);

	for (unsigned int i = 0; i < lines; i++)  
	{
		std::vector<Cell> init_lines;
		init_lines.clear();
		init_lines.reserve(cols);

		for (unsigned int j = 0; j < cols; j++)  
			init_lines.emplace_back(Cell());
		cells.emplace_back(init_lines);
	}
}

void init_locks(std::vector<std::vector<omp_lock_t>> &myLocks, long  lines, long cols)
{
    myLocks.clear();
	myLocks.reserve(lines);

    for (unsigned int i = 0; i < lines; i++)  
	{
		std::vector<omp_lock_t> init_lines(cols);
		
		for (unsigned int j = 0; j < cols; j++)
            omp_init_lock(&(init_lines[j]));
        
        myLocks.emplace_back(init_lines);
	}
}

std::vector<move> check_moves(item type, int i, int j)
{
    std::vector<move> moves;

    if (world[i][j].flag == 1)
        return moves;
    if (type == FOX)
    {
        if (i-1 >= 0)
            if (world_bck[i-1][j].resident.type == RABBIT)
                moves.insert(moves.end(), UP);
        if (j+1 < N)
            if (world_bck[i][j+1].resident.type == RABBIT)
                moves.insert(moves.end(), RIGHT);     
        if (i+1 < M)
            if (world_bck[i+1][j].resident.type == RABBIT)
                moves.insert(moves.end(), DOWN);      
        if (j-1 >= 0)
            if (world_bck[i][j-1].resident.type == RABBIT)
                moves.insert(moves.end(), LEFT);
    }
    if (type == RABBIT || type == FOX && moves.size() == 0)
    {
        if (i-1 >= 0)
            if (world_bck[i-1][j].resident.type == NONE)
                moves.insert(moves.end(), UP);
        if (j+1 < N)
            if (world_bck[i][j+1].resident.type == NONE)
                moves.insert(moves.end(), RIGHT);
        if (i+1 < M)
            if (world_bck[i+1][j].resident.type == NONE)
                moves.insert(moves.end(), DOWN);
        if (j-1 >= 0)
            if (world_bck[i][j-1].resident.type == NONE)
                moves.insert(moves.end(), LEFT);  
    }

    return moves;
}

std::pair<int, int> get_position(std::vector<move> moves, int k, int i, int j)
{
    int temp_i, temp_j;

    switch (moves[k])
    {
        case 0:
            temp_i = i-1;
            temp_j = j;
            break;
        case 1:
            temp_i = i;
            temp_j = j+1;
            break;
        case 2:
            temp_i = i+1;
            temp_j = j;
            break;
        case 3: 
            temp_i = i;
            temp_j = j-1;
            break;
        default:
            break;
    }   

    return std::make_pair(temp_i, temp_j);         
}

void print_world()
{
    printf("\n---------------------------------");
    printf("\n   ");
    for (int k=0; k<N; k++)
    {   
        if (k<10)
            printf("0%d|", k);
        else
            printf("%d|", k);
    }
    for (int i=0; i<M; i++)
    {
        printf("\n0%d:", i); 
        for (int j=0; j<N; j++)
        {          
            if (world[i][j].resident.type == FOX)
                printf(" F|");
            else if (world[i][j].resident.type == RABBIT)
                printf(" R|");
            else if (world[i][j].resident.type == ROCK)
                printf(" *|");
            else if (world[i][j].resident.type == NONE)
                printf("  |");
        } 
    }
    printf("\n---------------------------------");
}

void terminate()
{
    FOXES = RABBITS = ROCKS = 0;

    #pragma omp parallel for collapse(2) reduction (+ : FOXES, RABBITS, ROCKS)
    for (int i=0; i<M; i++)
    {    
        for(int j=0; j<N; j++)
        {    
            switch(world[i][j].resident.type)
            {
                case FOX:
                    FOXES++;
                    break;
                case RABBIT:
                    RABBITS++;
                    break;
                case ROCK:
                    ROCKS++;
                    break;
                default:
                    break;
            }
        }
    }
}

void kill_foxes(int target_age)
{        
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (world[i][j].resident.starving_age >= target_age && world[i][j].resident.type == FOX)
            {
                world[i][j].resident.type = NONE;
                world[i][j].resident.starving_age = 0;
                world[i][j].resident.breeding_age = 0;
            }
        }
    }
}

void destroy_locks(std::vector<std::vector<omp_lock_t>> &myLocks)
{
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            omp_destroy_lock(&(myLocks[i][j]));
        }
    }
}

void run_red(std::vector<std::vector<omp_lock_t>> &myLocks, int i, int j)
{   
    
    
    std::vector<move> moves = check_moves(world[i][j].resident.type, i, j);
    world[i][j].resident.breeding_age++;
    
    if (moves.size() > 0)
    {
        int new_move = (i * N + j) % moves.size();     
        std::pair <int, int> p = get_position(moves, new_move, i, j);

        int temp_i = p.first;
        int temp_j = p.second;

        // Check if the line being written is critical
        int set_lock = 0;

        int num_ts = omp_get_num_threads();
        int t_num = omp_get_thread_num();

        int rest =  M%num_ts;
        int extra = (t_num+1) <= rest ? 1 : 0;

        int small_pack = (int)floor(M/num_ts);      // used to calculate critical zone for divisions with rest
        int big_pack = small_pack + 1;

        int next_zone, previous_zone;
        if (extra){
            next_zone = (1+t_num)*big_pack;
            previous_zone = (t_num)*big_pack;
        }else{
            next_zone = (big_pack*rest)+(((1+t_num)-rest)*small_pack);
            previous_zone = (big_pack*rest)+(((t_num)-rest)*small_pack);
        }

        if( t_num!=(num_ts-1) && (temp_i == (next_zone-1) || temp_i == next_zone )){   // check conflict with next thread
            omp_set_lock(&(myLocks[temp_i][temp_j]));
            set_lock=1;
        }
        else if (t_num!=0 && (temp_i == (previous_zone-1) || temp_i == previous_zone )){  // check conflict with previous thread
            omp_set_lock(&(myLocks[temp_i][temp_j]));
            set_lock=1;
        }
            
        
        switch (world[i][j].resident.type)
        {
            case FOX:    
            
                //breeds?
                if (world[i][j].resident.breeding_age >= FOX_B_AGE)   
                {
                    world[i][j].resident.type = FOX;
                    world[i][j].resident.breeding_age = 0;
                    world[i][j].resident.starving_age = 0;
                }
                else
                    world[i][j].resident.type = NONE;
                
                // finds FOX in the destination cell -> FOX-FOX conflitc
                if (world[temp_i][temp_j].resident.type == FOX)
                {
                    if (world[i][j].resident.breeding_age > world[temp_i][temp_j].resident.breeding_age)
                    {
                        world[temp_i][temp_j].resident.breeding_age = world[i][j].resident.breeding_age;
                        world[temp_i][temp_j].resident.starving_age = world_bck[i][j].resident.starving_age + 1; 
                    }
                    else if (world[i][j].resident.breeding_age == world[temp_i][temp_j].resident.breeding_age)
                    { 
                        world[temp_i][temp_j].resident.starving_age = 
                            world_bck[i][j].resident.starving_age + 1 >= world[temp_i][temp_j].resident.starving_age ? 
                            world[temp_i][temp_j].resident.starving_age : world_bck[i][j].resident.starving_age + 1;
                    }                                 
                }
                // finds RABBIT or finds it EMPTY
                else 
                {                        
                    if (world[temp_i][temp_j].resident.type == RABBIT)
                        world[temp_i][temp_j].resident.starving_age = 0;
                    else
                    {
                        world[temp_i][temp_j].resident.starving_age = world_bck[i][j].resident.starving_age+1;
                    }

                    world[temp_i][temp_j].resident.type = FOX;
                    world[temp_i][temp_j].resident.breeding_age = world[i][j].resident.breeding_age;      
                }
            break;

            case RABBIT:
                // breeds?
                if (world[i][j].resident.breeding_age >= RABBIT_B_AGE)
                {
                    world[i][j].resident.type = RABBIT;
                    world[i][j].resident.breeding_age = 0;
                }
                else 
                    world[i][j].resident.type = NONE;
                // what is on the destination position?
                switch (world[temp_i][temp_j].resident.type)
                {
                    case FOX:
                        world[temp_i][temp_j].resident.starving_age = 0;
                        break;

                    case RABBIT:
                        if (world[i][j].resident.breeding_age > world[temp_i][temp_j].resident.breeding_age)
                            world[temp_i][temp_j].resident.breeding_age = world[i][j].resident.breeding_age;
                        break;

                    case NONE:
                    default:
                        world[temp_i][temp_j].resident.breeding_age = world[i][j].resident.breeding_age;     
                        world[temp_i][temp_j].resident.type = RABBIT;
                        break;                                                        
                }                                                     
            break;
        default:
            break;

        }

        world[temp_i][temp_j].flag = 1;
        moves.clear();

        if (set_lock){
            omp_unset_lock(&(myLocks[temp_i][temp_j]));
            set_lock=0;
        }

    }              
    
    else if (moves.size() == 0)
        ++world[i][j].resident.starving_age;
}

void run_black( std::vector<std::vector<omp_lock_t>> &myLocks, int i, int j)
{   
    std::vector<move> moves = check_moves(world[i][j].resident.type, i, j); 

    if (world[i][j].flag == 0)
        world[i][j].resident.breeding_age++;

    if (moves.size() > 0)
    {
        int new_move = (i * N + j) % moves.size();
        std::pair <int, int> p = get_position(moves, new_move, i, j);

        int temp_i = p.first;
        int temp_j = p.second;
        
        // Check if the line being written is critical
        int set_lock = 0;

        int num_ts = omp_get_num_threads();
        int t_num = omp_get_thread_num();

        int rest =  M%num_ts;
        int extra = (t_num+1) <= rest ? 1 : 0;

        int small_pack = (int)floor(M/num_ts); // used to calculate critical zone for divisions with rest
        int big_pack = small_pack + 1;

        int next_zone, previous_zone;
        if (extra){
            next_zone = (1+t_num)*big_pack;
            previous_zone = (t_num)*big_pack;
        }else{
            next_zone = (big_pack*rest)+(((1+t_num)-rest)*small_pack);
            previous_zone = (big_pack*rest)+(((t_num)-rest)*small_pack);
        }

        if( t_num!=(num_ts-1) && (temp_i == (next_zone-1) || temp_i == next_zone )){   // check conflict with next thread
            omp_set_lock(&(myLocks[temp_i][temp_j]));
            set_lock=1;
        }
        else if (t_num!=0 && (temp_i == (previous_zone-1) || temp_i == previous_zone )){  // check conflict with previous thread
            omp_set_lock(&(myLocks[temp_i][temp_j]));
            set_lock=1;
        }

        switch (world[i][j].resident.type)
        {
            case FOX:  

                // breeds?
                if (world[i][j].resident.breeding_age >= FOX_B_AGE)   
                {
                    world[i][j].resident.type = FOX;
                    world[i][j].resident.breeding_age = 0;
                    world[i][j].resident.starving_age = 0;
                }
                else
                    world[i][j].resident.type = NONE;

                // finds FOX in the destination cell -> FOX-FOX conflitc
                if (world[temp_i][temp_j].resident.type == FOX)
                {
                    if (world[i][j].resident.breeding_age > world[temp_i][temp_j].resident.breeding_age)
                    {
                        world[temp_i][temp_j].resident.breeding_age = world[i][j].resident.breeding_age;
                        world[temp_i][temp_j].resident.starving_age = world_bck[i][j].resident.starving_age + 1; 
                    }
                    else if (world[i][j].resident.breeding_age == world[temp_i][temp_j].resident.breeding_age)
                    { 
                        world[temp_i][temp_j].resident.starving_age = 
                            world_bck[i][j].resident.starving_age + 1 >= world[temp_i][temp_j].resident.starving_age ? 
                            world[temp_i][temp_j].resident.starving_age : world_bck[i][j].resident.starving_age + 1;
                    }                                            
                } 
                // finds RABBIT or finds it EMPTY
                else 
                {                            
                    if (world[temp_i][temp_j].resident.type == RABBIT)
                        world[temp_i][temp_j].resident.starving_age = 0;
                    else
                        world[temp_i][temp_j].resident.starving_age = world_bck[i][j].resident.starving_age+1;
                    
                    world[temp_i][temp_j].resident.type = FOX;
                    world[temp_i][temp_j].resident.breeding_age = world[i][j].resident.breeding_age;  
                }  
            break;

            case RABBIT:

                // breeds?
                if (world[i][j].resident.breeding_age >= RABBIT_B_AGE)
                {
                    world[i][j].resident.type = RABBIT;
                    world[i][j].resident.breeding_age = 0;
                }
                else 
                    world[i][j].resident.type = NONE;

                // what is on the destination position?
                switch (world[temp_i][temp_j].resident.type)
                {
                    case FOX:
                        world[temp_i][temp_j].resident.starving_age = 0;
                        break;
                    case RABBIT:
                        if (world[i][j].resident.breeding_age > world[temp_i][temp_j].resident.breeding_age)
                            world[temp_i][temp_j].resident.breeding_age = world[i][j].resident.breeding_age;
                        break;
                    case NONE:
                    default:
                        world[temp_i][temp_j].resident.type = RABBIT;
                        world[temp_i][temp_j].resident.breeding_age = world[i][j].resident.breeding_age; 
                        break;                                                        
                }         
            break;

        default:
            break;

        }
        moves.clear();

        if (set_lock){
            omp_unset_lock(&(myLocks[temp_i][temp_j]));
            set_lock=0;
        }
    }               
    
    else if (world[i][j].flag == 0 and moves.size() == 0) 
        ++world[i][j].resident.starving_age;

    else if (world[i][j].flag == 1)
        world[i][j].flag = 0;
}


void backup_world() { std::copy(std::begin(world), std::end(world), std::begin(world_bck)); }


int main(int argc, char *argv[])
{
    GENERATIONS = std::stol(argv[1]);
    M = std::stol(argv[2]);
    N = std::stol(argv[3]);

    N_ROCKS = std::stol(argv[4]);;
    N_RABBITS = std::stol(argv[5]);;
    N_FOXES = std::stol(argv[7]);;    

    RABBIT_B_AGE = std::stol(argv[6]);
    FOX_B_AGE = std::stol(argv[8]);
    FOX_S_AGE = std::stol(argv[9]);

    uint32_t seed = atoi(argv[10]);
    int cursor = 0;
    

    init_world(world, M, N);
    init_world(world_bck, M, N);

    std::vector<std::vector<omp_lock_t>> myLocks;
    init_locks(myLocks, M, N);

    generate_element(N_ROCKS, &seed, ROCK, 0, 0);
    generate_element(N_RABBITS, &seed, RABBIT, 0, 0);
    generate_element(N_FOXES, &seed, FOX, 0, 0);

    backup_world();

    
    double exec_time;
    exec_time = -omp_get_wtime();
    
    while (cursor < GENERATIONS)
    {   
        // red gen
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < M; i++){
            for (int j = (i % 2 == 0) ? 0 : 1; j < N; j = j+2){
                run_red(myLocks, i, j);
            }
        }
 
        backup_world();  

        // black gen
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < M; i++){
            for (int j = (i % 2 == 0) ? 1 : 0; j < N; j = j+2){
                run_black(myLocks, i, j);
            }
        }
        
        kill_foxes(FOX_S_AGE);                          
        backup_world();   
        cursor++;
    }

    terminate();
    
    exec_time += omp_get_wtime();
    destroy_locks(myLocks);

    fprintf(stderr, "%.1fs", exec_time); fflush(stdout);
    printf("\n%d %d %d\n", ROCKS, RABBITS, FOXES); fflush(stdout);

    return 0;
}
