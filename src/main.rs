use std::env;
use numext_fixed_uint::U512;
use rayon::prelude::IntoParallelIterator;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use boolean_function_extender::BooleanFunctionTester;
use boolean_function_extender::u32_tester::U32Tester;
use boolean_function_extender::u512_tester::U512Tester;

const RING_SIZE: usize = 9;
const RULE_NUMBER: u32 = 1438886595;

fn main() {
    /*let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: {} <rounds_count>", args[0]);
        return;
    }
    let rounds_count: usize = match args[1].parse() {
        Ok(value) => value,
        Err(_) => {
            eprintln!("Invalid rounds count, must be a positive integer");
            return;
        }
    };

    let mut total_res_sac = 0;

    for round_number in 0..rounds_count {
        println!("Round {}", round_number);

    }*/

    println!("Generating extended functions");
    let extended_booleans: Vec<(u32, U512)> = (0..=u32::MAX/100).into_par_iter().map(|rule_number| {
        (rule_number, extend_rule_5_to_9(rule_number))
    }).collect();
    println!("Generated extended functions");

    println!("Testing strict avalanche criterion (SAC)");
    let res_sac: usize = extended_booleans.par_iter().map(|(rule_number, output_9_rule_number)| {
        if U512Tester::is_strict_avalanche_criterion_ok(output_9_rule_number) == U32Tester::is_strict_avalanche_criterion_ok(rule_number) {
            1
        } else {
            0
        }
    }).sum();
    println!("equal {}", res_sac);

    println!("Testing first order correlation immunity");
    let res_foci: usize = extended_booleans.par_iter().map(|(rule_number, output_9_rule_number)| {
        if U512Tester::is_first_order_correlation_immune(output_9_rule_number) == U32Tester::is_first_order_correlation_immune(rule_number) {
            1
        } else {
            0
        }
    }).sum();
    println!("equal {}", res_foci);

    println!("Testing balanced");
    let res_balanced: usize = extended_booleans.par_iter().map(|(rule_number, output_9_rule_number)| {
        if U512Tester::is_function_balanced(output_9_rule_number) == U32Tester::is_function_balanced(rule_number) {
            1
        } else {
            0
        }
    }).sum();
    println!("equal {}", res_balanced);

    println!("Testing propagation criterion deg 2");
    let res_propagation_2: usize = extended_booleans.par_iter().map(|(rule_number, output_9_rule_number)| {
        if U512Tester::is_propagation_criterion_deg_k_ok(output_9_rule_number, 2) == U32Tester::is_propagation_criterion_deg_k_ok(rule_number, 2) {
            1
        } else {
            0
        }
    }).sum();
    println!("equal {}", res_propagation_2);

    println!("Testing propagation criterion deg 3");
    let res_propagation_3: usize = extended_booleans.par_iter().map(|(rule_number, output_9_rule_number)| {
        if U512Tester::is_propagation_criterion_deg_k_ok(output_9_rule_number, 3) == U32Tester::is_propagation_criterion_deg_k_ok(rule_number, 3) {
            1
        } else {
            0
        }
    }).sum();
    println!("equal {}", res_propagation_3);

    println!("Testing propagation criterion deg 4");
    let res_propagation_4: usize = extended_booleans.par_iter().map(|(rule_number, output_9_rule_number)| {
        if U512Tester::is_propagation_criterion_deg_k_ok(output_9_rule_number, 4) == U32Tester::is_propagation_criterion_deg_k_ok(rule_number, 4) {
            1
        } else {
            0
        }
    }).sum();
    println!("equal {}", res_propagation_4);

    println!("Testing algebraic degree");
    let res_degree: usize = extended_booleans.par_iter().map(|(rule_number, output_9_rule_number)| {
        if U512Tester::get_function_degree(output_9_rule_number) == U32Tester::get_function_degree(rule_number) {
            1
        } else {
            0
        }
    }).sum();
    println!("equal {}", res_degree);
}

fn extend_rule_5_to_9(rule_number: u32) -> U512 {
    let mut output_rule_number = U512::zero();
    for i in 0usize..(1 << RING_SIZE) {
        let input_bits = unsigned_to_bool_array::<RING_SIZE>(i);
        let round_1 = get_new_ring(input_bits, rule_number);
        let round_2 = get_new_ring(round_1, rule_number);
        if round_2[4] {
            output_rule_number |= U512::one() << i;
            //print!("1");
        } else {
            //print!("0");
        }
    }
    output_rule_number
}



#[inline(always)]
fn compute_ca_rule(rule_number: u32, input_bits: u8) -> bool {
    return (rule_number & (1 << input_bits)) != 0;
}

fn get_new_ring(ring: [bool; RING_SIZE], rule_number: u32) -> [bool; RING_SIZE] {
    let mut new_ring = [false; RING_SIZE];
    for i in 0..RING_SIZE {
        let input_bits = (ring[(i + RING_SIZE - 2) % RING_SIZE] as u8) << 4
            | (ring[(i + RING_SIZE - 1) % RING_SIZE] as u8) << 3
            | (ring[i] as u8) << 2
            | (ring[(i + 1) % RING_SIZE] as u8) << 1
            | ring[(i + 2) % RING_SIZE] as u8;
        new_ring[i] = compute_ca_rule(rule_number, input_bits);
    }
    new_ring
}

#[inline(always)]
fn unsigned_to_bool_array<const S: usize>(number: usize) -> [bool; S] {
    let mut bits = [false; S];
    for i in 0..S {
        bits[i] = (number & (1 << i)) != 0;
    }
    bits
}

#[inline(always)]
fn bool_array_to_unsigned<const S: usize>(bits: [bool; S]) -> usize {
    let mut number = 0;
    for i in 0..S {
        if bits[i] {
            number |= 1 << i;
        }
    }
    number
}