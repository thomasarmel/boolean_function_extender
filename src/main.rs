use numext_fixed_uint::U512;
use rayon::prelude::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use boolean_function_extender::BooleanFunctionTester;
use boolean_function_extender::u32_tester::U32Tester;
use boolean_function_extender::u512_tester::U512Tester;

const RING_SIZE: usize = 9;
const RULE_NUMBER: u32 = 100/*1438886595*/;

fn main() {
    let res = (0../*=u32::MAX*/10000000u32).into_par_iter().map(|rule_number| {
        let output_9_rule_number = extend_rule_5_to_9(rule_number);
        if U32Tester::is_propagation_criterion_deg_k_ok(rule_number, 3) {
            println!("{} {}", rule_number, output_9_rule_number);
        }
        if U512Tester::is_propagation_criterion_deg_k_ok(output_9_rule_number, 2) == U32Tester::is_propagation_criterion_deg_k_ok(rule_number, 2) {
            1
        } else {
            0
        }
    }).count();
    println!("equal {}", res);
    let output_rule_number = extend_rule_5_to_9(RULE_NUMBER);
    println!("{}", output_rule_number);
    let anf = U512Tester::fast_bool_anf_transform_unsigned(output_rule_number.clone(), 9);
    println!("{}", anf);
    println!("degree {}", U512Tester::get_function_degree(output_rule_number));
    println!("original {}", U32Tester::get_function_degree(RULE_NUMBER));
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