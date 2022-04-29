#include <vector>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <unistd.h>


enum item {FOX, RABBIT, ROCK, NONE};
enum move {UP, RIGHT, DOWN, LEFT};

struct Resident
{
    item type;
    int starving_age;
    int breeding_age;

    Resident() : type(NONE), starving_age(0), breeding_age(0) {}
};

struct Conflict
{
    int start_i;
    int start_j;
    int dest_i;
    int dest_j;

    Conflict(int i,int j, int k,int l) {
        this->start_i = i;
        this->start_j = j;
        this->dest_i = k;
        this->dest_j = l;

    }
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
int LINE_BYTES;

int base_chunk, big_chunk, small_chunk;
int rest, p;

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
        printf("%02d|", k);
    }
    for (int i=0; i<M; i++)
    {
        printf("\n%02d:", i); 
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

void run_red(std::vector<std::vector<omp_lock_t>> &myLocks, int i, int j, int r)
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

        int rest =  r%num_ts;
        int extra = (t_num+1) <= rest ? 1 : 0;

        int small_pack = (int)floor(r/num_ts);      // used to calculate critical zone for divisions with rest
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

void run_black(std::vector<std::vector<omp_lock_t>> &myLocks, int i, int j, int r)
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

        int rest =  r%num_ts;
        int extra = (t_num+1) <= rest ? 1 : 0;

        int small_pack = (int)floor(r/num_ts); // used to calculate critical zone for divisions with rest
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


void comunication(MPI_Comm newworld, int p, int id, int s_i, int f_i, int * responsible) {


    std::vector<Cell> first_send_buffer_prev;
    std::vector<Cell> first_send_buffer_next;
    std::vector<Cell> first_rcv_buffer_next;
    std::vector<Cell> first_rcv_buffer_prev;
    std::vector<Cell> second_send_buffer_prev;
    std::vector<Cell> second_send_buffer_next;
    std::vector<Cell> second_rcv_buffer_next;
    std::vector<Cell> second_rcv_buffer_prev;

    // Individual buffers because of individual requests
    first_send_buffer_prev.reserve(N);
    first_send_buffer_next.reserve(N);
    first_rcv_buffer_next.reserve(N);
    first_rcv_buffer_prev.reserve(N);
    second_send_buffer_prev.reserve(N);
    second_send_buffer_next.reserve(N);
    second_rcv_buffer_next.reserve(N);
    second_rcv_buffer_prev.reserve(N);
    
    
    MPI_Request r_first_send_next, r_first_send_prev, r_first_rcv_next, r_first_rcv_prev, r_second_send_next, r_second_send_prev, r_second_rcv_next, r_second_rcv_prev;
    MPI_Status s_first_send_next, s_first_send_prev, s_first_rcv_next, s_first_rcv_prev, s_second_send_next, s_second_send_prev, s_second_rcv_next, s_second_rcv_prev;
    int first_line_under = s_i-1;
    int second_line_under = s_i-2;
    int first_line_over = f_i+1;
    int second_line_over = f_i+2;


    if(first_line_under>=0){
        
        for(int j=0; j<N; j++){
            first_send_buffer_prev[j]=world[s_i][j];
        }
        MPI_Irecv(&(first_rcv_buffer_prev[0]), LINE_BYTES, MPI_BYTE, responsible[first_line_under], id, newworld, &r_first_rcv_prev);
        MPI_Isend(&(first_send_buffer_prev[0]), LINE_BYTES, MPI_BYTE, responsible[first_line_under], responsible[first_line_under], newworld, &r_first_send_prev);

        MPI_Wait(&r_first_rcv_prev, &s_first_rcv_prev);
        //printf("\nid:%d, receiving line:%d\n", id, s_i-1);
        for(int j=0; j<N; j++){
            world[s_i-1][j]=first_rcv_buffer_prev[j];
        }
        if((s_i + 1)<= f_i){
            for(int j=0; j<N; j++){
                second_send_buffer_prev[j]=world[s_i+1][j];
            }
            MPI_Isend(&(second_send_buffer_prev[0]), LINE_BYTES, MPI_BYTE, responsible[first_line_under], responsible[first_line_under], newworld, &r_second_send_prev);
        }
        if(second_line_under>=0){
            MPI_Irecv(&(second_rcv_buffer_prev[0]), LINE_BYTES, MPI_BYTE, responsible[second_line_under], id, newworld, &r_second_rcv_prev);
            
            if (responsible[second_line_under] != responsible[first_line_under]){
                MPI_Wait(&r_first_send_prev, &s_first_send_prev);
                MPI_Isend(&(first_send_buffer_prev[0]), LINE_BYTES, MPI_BYTE, responsible[second_line_under], responsible[second_line_under], newworld, &r_first_send_prev);
            }

            MPI_Wait(&r_second_rcv_prev, &s_second_rcv_prev);
            for(int j=0; j<N; j++){
                world[s_i-2][j]=second_rcv_buffer_prev[j];
            }
            

        }
        
        
    }
    if(first_line_over<M){

        for(int j=0; j<N; j++){
            first_send_buffer_next[j]=world[f_i][j];
        }
        
        MPI_Irecv(&(first_rcv_buffer_next[0]), LINE_BYTES, MPI_BYTE, responsible[first_line_over], id, newworld, &r_first_rcv_next);
        MPI_Isend(&(first_send_buffer_next[0]), LINE_BYTES, MPI_BYTE, responsible[first_line_over], responsible[first_line_over], newworld, &r_first_send_next);
        
        MPI_Wait(&r_first_rcv_next, &s_first_rcv_next);
        
        for(int j=0; j<N; j++){
            world[f_i+1][j]=first_rcv_buffer_next[j];
        }
        if((f_i - 1)>= s_i){
            for(int j=0; j<N; j++){
                second_send_buffer_next[j]=world[f_i-1][j];
            }
            MPI_Isend(&(second_send_buffer_next[0]), LINE_BYTES, MPI_BYTE, responsible[first_line_over], responsible[first_line_over], newworld, &r_second_send_next);
        }
        if(second_line_over<M){
            MPI_Irecv(&(second_rcv_buffer_next[0]), LINE_BYTES, MPI_BYTE, responsible[second_line_over], id, newworld, &r_second_rcv_next);
            
            if (responsible[second_line_over] != responsible[first_line_over]){
                MPI_Wait(&r_first_send_next, &s_first_send_next);
                MPI_Isend(&(first_send_buffer_next[0]), LINE_BYTES, MPI_BYTE, responsible[second_line_over], responsible[second_line_over], newworld, &r_first_send_next);
            }

            MPI_Wait(&r_second_rcv_next, &s_second_rcv_next);
            for(int j=0; j<N; j++){
                world[f_i+2][j]=second_rcv_buffer_next[j];
            }
            

        }
        
        
    }

    

}

int main(int argc, char *argv[])
{   
    MPI_Init (&argc, &argv);

    GENERATIONS = std::stol(argv[1]);
    M = std::stol(argv[2]);
    N = std::stol(argv[3]);
    LINE_BYTES = (sizeof(Cell) * N);

    N_ROCKS = std::stol(argv[4]);;
    N_RABBITS = std::stol(argv[5]);;
    N_FOXES = std::stol(argv[7]);;    

    RABBIT_B_AGE = std::stol(argv[6]);
    FOX_B_AGE = std::stol(argv[8]);
    FOX_S_AGE = std::stol(argv[9]);

    uint32_t seed = atoi(argv[10]);
    int cursor = 0;

    MPI_Status status;
    int id, p, chunk;
    int s_i, f_i, r;

    
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    MPI_Comm newworld;
    MPI_Group new_group;
    MPI_Comm_group(MPI_COMM_WORLD, &new_group);
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &newworld);

    int process_limit = p;

    // remove too much processes
    while( (p>1) && (p>M) ){

        process_limit--;
        // Remove all unnecessary ranks
        int ranges[1][3] = {{ process_limit, p-1, 1 }};
        MPI_Group_range_excl(new_group, 1, ranges, &new_group);


        MPI_Comm_create(newworld, new_group, &newworld);

        if (newworld == MPI_COMM_NULL)
        {
        MPI_Finalize();
        exit(0);
        }

        MPI_Comm_size (newworld, &p);
    }

    base_chunk = (int)floor(M/p);    
    rest = M%p;

    big_chunk = base_chunk+1;
    small_chunk = base_chunk;

    int extra, start, finish, range;
    int *responsible =(int*) malloc(M*sizeof(int));


    init_world(world, M, N);
    init_world(world_bck, M, N);

    std::vector<std::vector<omp_lock_t>> myLocks;
    init_locks(myLocks, M, N);

    generate_element(N_ROCKS, &seed, ROCK, 0, 0);
    generate_element(N_RABBITS, &seed, RABBIT, 0, 0);
    generate_element(N_FOXES, &seed, FOX, 0, 0);

    backup_world();
    double exec_time;
    if(!id){
        exec_time = -MPI_Wtime();
    }

    extra = (id+1) <= rest ? 1 : 0;
    chunk = base_chunk + extra;

    int temp_rest = rest;
    int temp_id = 0;
    int temp_chunk = 0;
    int line = 0;
    while(line<M){
        temp_chunk= temp_rest > 0 ? big_chunk : small_chunk;
        while(temp_chunk>0){
            responsible[line] = temp_id;
            temp_chunk--;
            line++;
        }
        temp_rest--;
        temp_id++;
    }

    if (extra){
        s_i = (id)*big_chunk;
    }else{
        s_i = (big_chunk*rest)+(((id)-rest)*small_chunk);
    }

    f_i = s_i + (chunk -1);

    start=std::max(0, s_i-1);
    finish=std::min(M-1, f_i+1);
    range = (finish - start) + 1;

    
    while (cursor < GENERATIONS)
    {   
        /*
        ****************************************************************************************************************************************
        ***********************************************************                *************************************************************
        **********************************************************       RED        ***********************************************************
        ***********************************************************                *************************************************************
        ****************************************************************************************************************************************
        */
        // red gen
        #pragma omp parallel for schedule(static)
        for (int i = start; i <= finish; i++){
            for (int j = (i % 2 == 0) ? 0 : 1; j < N; j = j+2){
                run_red(myLocks, i, j, range);
            }
        }
        comunication(newworld, p, id, s_i, f_i, responsible);
        /*
        ****************************************************************************************************************************************
        ***********************************************************                *************************************************************
        **********************************************************      BLACK       ***********************************************************
        ***********************************************************                *************************************************************
        ****************************************************************************************************************************************
        */
        backup_world(); 

        // black gen
        #pragma omp parallel for schedule(static)
        for (int i = start; i <= finish; i++){
            for (int j = (i % 2 == 0) ? 1 : 0; j < N; j = j+2){
                run_black(myLocks, i, j, range);
            }
        }
        comunication(newworld, p, id, s_i, f_i, responsible);

        kill_foxes(FOX_S_AGE);                         
        backup_world();   
        cursor++;
    }

    /*
    ****************************************************************************************************************************************
    ***********************************************************                *************************************************************
    *********************************************************    FINALIZATION    ***********************************************************
    ***********************************************************                *************************************************************
    ****************************************************************************************************************************************
    */

    std::vector<Cell> send_buffer;
    send_buffer.reserve(N*chunk);

    int scount = chunk * LINE_BYTES;

    int i=0;
    for(int k = s_i; k<=f_i; k++){
        for(int j=0; j<N; j++){
            send_buffer[i]=world[k][j];
            i++;
        }
    }

    if(!id){
        std::vector<Cell> rcv_buffer;
        int * chunks = (int*) malloc(p*sizeof(int));

        /* displacement vectors for communication */
        int *rdispls=(int*) malloc(p*sizeof(int));
        int *rdisplsB=(int*) malloc(p*sizeof(int));
        int *rcounts=(int*) malloc(p*sizeof(int));

        
        int tmp_rest = rest;
        int displ = 0;
        for (int i=0; i<p;i++){

            chunks[i]= tmp_rest > 0 ? big_chunk : small_chunk;
            tmp_rest--;

            rdisplsB[i]=displ*LINE_BYTES;
            rdispls[i]=displ*N;
            rcounts[i]=(chunks[i])*LINE_BYTES;

            displ += (chunks[i]);
        }

        rcv_buffer.reserve(M*N);

        MPI_Gatherv(&(send_buffer[0]), scount, MPI_BYTE, &(rcv_buffer[0]), rcounts, rdisplsB, MPI_BYTE, 0, newworld);

        /* Unpacking of the info from the buffer */
        for(int i=0; i<p; i++){

            int chunk = chunks[i];
            int s_i, f_i;
            int extra = (i+1) <= rest ? 1 : 0;

            if (extra){
                s_i = (i)*big_chunk;
            }else{
                s_i = (big_chunk*rest)+(((i)-rest)*small_chunk);
            }

            f_i = s_i + (chunk -1);

            int inicio = rdispls[i];
            int j=0;
            for(int k=s_i; k<=f_i; k++){
                for(int l=0; l<N; l++){
                    world[k][l]=rcv_buffer[inicio + j];
                    j++;
                }
            }   
        }

        free(chunks);
        free(rdispls);
        free(rdisplsB);
        free(rcounts);
        terminate();
        exec_time += MPI_Wtime();
        fprintf(stderr, "%.1fs", exec_time); fflush(stdout);
        printf("\n%d %d %d\n", ROCKS, RABBITS, FOXES); fflush(stdout);
    }else{
        MPI_Gatherv(&(send_buffer[0]), scount, MPI_BYTE, NULL, NULL, NULL, MPI_BYTE, 0, newworld);
    }


    destroy_locks(myLocks);
    MPI_Finalize();
    return 0;
}