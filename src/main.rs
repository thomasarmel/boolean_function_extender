use fast_boolean_anf_transform::fast_bool_anf_transform_unsigned;
use numext_fixed_uint::U512;
use rayon::prelude::IntoParallelIterator;
use rayon::iter::ParallelIterator;
use boolean_function_extender::BooleanFunctionTester;
use boolean_function_extender::u32_tester::U32Tester;
use boolean_function_extender::u512_tester::U512Tester;

const RING_SIZE: usize = 9;
const EQUIVALENCE_CLASSES: [u32; 48] = [0xaa55aa55, 0xaa55ab55, 0xaa55bb55, 0xaa5dbb55, 0xaaddbb55, 0xaa5dbb51, 0x2a5dbb51, 0xaaddbb51, 0x2a5dbf51, 0x6a5dbb51, 0x2addbb51, 0xa8ddbb51, 0xaeddda51, 0x0a5dbf51, 0x8addda51, 0xa8dd9b51, 0x88ddbb51, 0x88ddbb11, 0x8c5dda51, 0xa89d9b51, 0x8eddda51, 0xaefdda51, 0x025dbf51, 0x88ddda51, 0x88dd9b51, 0xceddda51, 0x0eddda51, 0x425dbf51, 0x8cddda51, 0x88dddb51, 0x289d9b51, 0x86fdda51, 0x88dddb71, 0xcefdda51, 0x0efdda51, 0x288d9b51, 0x8cfdda51, 0x8cdddb51, 0x8ccdda51, 0x289d9b41, 0x488ddb51, 0xccfdda51, 0x688d9b51, 0x288d9b41, 0x288d1b41, 0xdcfdda51, 0x68ad9b51, 0x688ddb51];

fn main() {
    /*(0..=u32::MAX).into_par_iter().for_each(|rule_number| {
        let output_9_rule_number = extend_rule_5_to_9(rule_number);
        let initial_deg = U32Tester::get_function_degree(&rule_number);
        if U512Tester::get_function_degree(&output_9_rule_number) == initial_deg {
            println!("{} -> {} (degree = {})", rule_number, output_9_rule_number, initial_deg);
        }
    });*/
    (0..=u32::MAX).into_par_iter().for_each(|f|  {
        let walsh_spectrum = U32Tester::absolute_walsh_spectrum(&f);
        let autocorrelation_spectrum = U32Tester::absolute_autocorrelation_spectrum(&f);
        let equivalent_classes_count: u32 = EQUIVALENCE_CLASSES.iter().map(|eq| {
            let walsh_spectrum_eq = U32Tester::absolute_walsh_spectrum(&eq);
            let autocorrelation_spectrum_eq = U32Tester::absolute_autocorrelation_spectrum(&eq);
            if walsh_spectrum == walsh_spectrum_eq && autocorrelation_spectrum == autocorrelation_spectrum_eq {
                1
            } else {
                0
            }
        }).sum();
        if equivalent_classes_count != 1 {
            println!("{} -> {}", f, equivalent_classes_count);
        }
    });
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

#[allow(dead_code)]
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