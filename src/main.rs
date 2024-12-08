use rand::{random, Rng};
#[derive(Clone)]
struct NBoard {
    n: usize,
    queens: Vec<usize>,
}

impl NBoard {
    fn new_random(n: usize) -> Self {
        let queens: Vec<usize> = (0..n).map(|_| random::<usize>() % n).collect();
        NBoard { n, queens }
    }
    fn calculate_energy(&self) -> usize {
        let mut energy = 0;

        let mut seen_rows: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut seen_main_diagonals: std::collections::HashSet<isize> =
            std::collections::HashSet::new();
        let mut seen_anti_diagonals: std::collections::HashSet<isize> =
            std::collections::HashSet::new();

        for col in 0..self.n {
            let row = self.queens[col];
            if !seen_rows.insert(row) {
                energy += 1;
            }
            if !seen_main_diagonals.insert(row as isize - col as isize) {
                energy += 1;
            }
            if !seen_anti_diagonals.insert(row as isize + col as isize) {
                energy += 1;
            }
        }
        energy
    }
    fn perturb(&mut self) {
        let board_size = self.n;

        // randomly select a column
        let mut rng = rand::thread_rng();
        let col = rng.gen_range(0..board_size);

        // randomly selec a row for the queen
        let current_row = self.queens[col];
        let mut new_row: usize;
        loop {
            new_row = rng.gen_range(0..board_size);
            if new_row != current_row {
                break;
            }
        }
        self.queens[col] = new_row;
    }
    fn visualize(&self) {
        let n = self.queens.len();

        for row in 0..n {
            for col in 0..n {
                if self.queens[col] == row {
                    print!("Q ");
                } else {
                    print!(". ");
                }
            }
            println!();
        }
    }
    fn get_neighbors(&self) -> Vec<NBoard> {
        let mut neighbors = Vec::new();
        for col in 0..self.n {
            let mut new_board = self.clone();
            for row in 0..self.n {
                if new_board.queens[col] != row {
                    new_board.queens[col] = row;
                    neighbors.push(new_board.clone());
                }
            }
        }
        neighbors
    }
}
fn selection(population: &Vec<NBoard>, n: usize) -> NBoard {
    let mut rng = rand::thread_rng();
    let mut best = &population[0];
    let mut best_energy = best.calculate_energy();
    for _ in 0..n {
        let candidate = &population[rng.gen_range(0..population.len())];
        let candidate_energy = candidate.calculate_energy();
        if candidate_energy < best_energy {
            best = candidate;
            best_energy = candidate_energy;
        }
    }
    best.clone()
}
fn mutate(board: &mut NBoard) {
    board.perturb();
}
fn crossover(parent1: &NBoard, parent2: &NBoard) -> NBoard {
    let mut rng = rand::thread_rng();
    let crossover_point = rng.gen_range(0..parent1.n);

    let mut child_queens = parent1.queens[0..crossover_point].to_vec();
    child_queens.extend_from_slice(&parent2.queens[crossover_point..]);

    NBoard {
        n: parent1.n,
        queens: child_queens,
    }
}
fn genetic_algorithm(n: usize, population_size: usize, generations: usize) -> NBoard {
    let mut population: Vec<NBoard> = (0..population_size)
        .map(|_| NBoard::new_random(n))
        .collect();

    for generation in 0..generations {
        let mut new_population: Vec<NBoard> = Vec::new();

        // Selection and crossover
        while new_population.len() < population_size {
            let parent1 = selection(&population, 3); // Tournament selection
            let parent2 = selection(&population, 3);
            let child = crossover(&parent1, &parent2);
            new_population.push(child);
        }

        // Mutation
        for individual in new_population.iter_mut() {
            if rand::random::<f64>() < 0.1 {
                // 10% mutation rate
                mutate(individual);
            }
        }

        // Replace population with the new one
        population = new_population;

        // Check if we have a solution
        let best = population
            .iter()
            .min_by_key(|b| b.calculate_energy())
            .unwrap();
        if best.calculate_energy() == 0 {
            return best.clone();
        }
    }

    // Return the best solution after all generations
    population
        .iter()
        .min_by_key(|b| b.calculate_energy())
        .unwrap()
        .clone()
}

// implement simulated annealing algorithm
fn simulated_annealing(
    board: NBoard,
    initial_temp: f64,
    cooling_rate: f64,
    min_temp: f64,
) -> NBoard {
    let mut current_board = board.clone();
    let mut current_energy = current_board.calculate_energy();
    let mut temperature = initial_temp;

    while temperature > min_temp {
        let mut neighbor = current_board.clone();
        neighbor.perturb();
        let neighbor_energy = neighbor.calculate_energy();

        let energy_difference = neighbor_energy as isize - current_energy as isize;
        let acceptance_prob: f64 = (-energy_difference as f64 / temperature).exp();

        let mut rng = rand::thread_rng();
        if rng.gen::<f64>() < acceptance_prob {
            current_board = neighbor;
            current_energy = neighbor_energy;
        }
        temperature *= cooling_rate;
    }
    current_board
}
// hillclimbing
fn hill_climbing(board: NBoard) -> NBoard {
    let mut current_board = board;
    let mut current_energy = current_board.calculate_energy();

    loop {
        let neighbors = current_board.get_neighbors();
        if neighbors.is_empty() {
            break; // No more neighbors, we're done
        }

        let mut best_neighbor = neighbors[0].clone();
        let mut best_energy = best_neighbor.calculate_energy();

        for neighbor in neighbors {
            let energy = neighbor.calculate_energy();
            if energy < best_energy {
                best_energy = energy;
                best_neighbor = neighbor;
            }
        }
        if best_energy >= current_energy {
            break;
        }

        current_board = best_neighbor;
        current_energy = best_energy;
    }

    current_board
}

//helper
fn display_solution(board: NBoard, problem_size: usize) {
    println!("Solution for {}-queen problem:", problem_size);
    board.visualize();
    println!("Final energy: {}", board.calculate_energy());
}
fn main() {
    println!("\n==========S ANNEAHLING===========\n");

    /// display anneahling  
    let initial_temp = 10000.0;
    let cooling_rate = 0.95;
    let min_temp = 0.01;

    let n_4 = 4;
    let board_4 = NBoard::new_random(n_4);
    println!("Initial board for {}-queen problem:", n_4);
    board_4.visualize();
    println!("Initial energy: {}", board_4.calculate_energy());

    let solution_4 = simulated_annealing(board_4, initial_temp, cooling_rate, min_temp);
    display_solution(solution_4, n_4);

    println!("\n=====================\n");

    let n_8 = 8;
    let board_8 = NBoard::new_random(n_8);
    println!("Initial board for {}-queen problem:", n_8);
    board_8.visualize();
    println!("Initial energy: {}", board_8.calculate_energy());

    let solution_8 = simulated_annealing(board_8, initial_temp, cooling_rate, min_temp);
    display_solution(solution_8, n_8);
    println!("\n==========HILL CLIMBING===========\n");

    // display hill climbing
    let n_4 = 4;
    let board_4 = NBoard::new_random(n_4);
    println!("Initial board for {}-queen problem:", n_4);
    board_4.visualize();
    println!("Initial energy: {}", board_4.calculate_energy());

    let solution_4 = hill_climbing(board_4);
    display_solution(solution_4, n_4);

    println!("\n=====================\n");

    let n_8 = 8;
    let board_8 = NBoard::new_random(n_8);
    println!("Initial board for {}-queen problem:", n_8);
    board_8.visualize();
    println!("Initial energy: {}", board_8.calculate_energy());

    let solution_8 = hill_climbing(board_8);
    display_solution(solution_8, n_8);
    println!("\n=========GENETIC============\n");

    // genetic algorithm
    let n_4 = 4;
    let n_8 = 8;
    let population_size = 100;
    let generations = 1000;

    let solution_4 = genetic_algorithm(n_4, population_size, generations);
    display_solution(solution_4, n_4);

    println!("\n=====================\n");

    let solution_8 = genetic_algorithm(n_8, population_size, generations);
    display_solution(solution_8, n_8);
}
