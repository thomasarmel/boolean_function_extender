use std::arch::x86_64::{_bzhi_u32, _popcnt32};
use gen_combinations::CombinationIterator;
use crate::BooleanFunctionTester;

pub struct U32Tester;

impl BooleanFunctionTester for U32Tester {
    type UnsignedRepr = u32;
    const NUM_VARIABLES: usize = 5;

    fn fast_bool_anf_transform_unsigned(rule_number: &Self::UnsignedRepr, num_variables_function: usize) -> Self::UnsignedRepr {
        fast_boolean_anf_transform::fast_bool_anf_transform_unsigned(*rule_number, num_variables_function)
    }

    fn get_function_degree(rule_number: &Self::UnsignedRepr) -> usize {
        let anf_form = Self::fast_bool_anf_transform_unsigned(rule_number, Self::NUM_VARIABLES);
        let max_input_value = unsafe { _bzhi_u32(u32::MAX, Self::NUM_VARIABLES as u32) };
        (0..=max_input_value).into_iter().map(|bit_position| {
            if anf_form & (1 << bit_position) != 0 {
                bit_position.count_ones() as usize
            } else {
                0
            }
        }).max().unwrap_or(0)
    }

    fn is_strict_avalanche_criterion_ok(rule_number: &Self::UnsignedRepr) -> bool {
        let max_input_value = unsafe { _bzhi_u32(u32::MAX, Self::NUM_VARIABLES as u32) };
        (0..Self::NUM_VARIABLES).into_iter().all(|constant_position| {
            let constant = 1 << constant_position;
            (0..=max_input_value).into_iter().filter(|&x| {
                let x_prime = x ^ constant;
                Self::compute_cellular_automata_rule(rule_number, x) == Self::compute_cellular_automata_rule(rule_number, x_prime)
            }).count() == (1 << (Self::NUM_VARIABLES - 1))
        })
    }

    fn compute_cellular_automata_rule(rule_number: &Self::UnsignedRepr, input_bits: u32) -> bool {
        let max_input_value = unsafe { _bzhi_u32(u32::MAX, Self::NUM_VARIABLES as u32) };
        #[cfg(debug_assertions)]
        if input_bits > max_input_value {
            panic!("Input bits must be less or equal than {}", max_input_value);
        }
        (rule_number & (1 << input_bits)) != 0
    }

    fn is_function_balanced(rule_number: &Self::UnsignedRepr) -> bool {
        const EXPECTED_SET_NUMBER: i32 = 1 << 4;
        unsafe {
            _popcnt32(*rule_number as i32) == EXPECTED_SET_NUMBER
        }
    }

    fn is_first_order_correlation_immune(rule_number: &Self::UnsignedRepr) -> bool {
        (0..Self::NUM_VARIABLES)
            .map(|input_bit_number| {
                1 << input_bit_number
            })
            .all(|w| {
                Self::fast_walsh_transform(rule_number, w) == 0
            })
    }

    fn fast_walsh_transform(rule_number: &Self::UnsignedRepr, w: u32) -> i32 {
        let max_input_value = unsafe { _bzhi_u32(u32::MAX, Self::NUM_VARIABLES as u32) };
        (0..=max_input_value).map(|x| {
            if Self::compute_cellular_automata_rule(rule_number, x) {
                if Self::fast_binary_dot_product(w, x as u32) & 1 == 0 { // % modulo 2
                    1
                } else {
                    -1
                }
            } else {
                0
            }
        }).sum()
    }

    fn is_propagation_criterion_deg_k_ok(rule_number: &Self::UnsignedRepr, k: usize) -> bool {
        let max_input_value = unsafe { _bzhi_u32(u32::MAX, Self::NUM_VARIABLES as u32) };
        if k == 0 {
            return true;
        }
        let strict_avalanche_criterion_ok = Self::is_strict_avalanche_criterion_ok(rule_number);
        if k == 1 {
            return strict_avalanche_criterion_ok;
        }
        if k >= 1 && !strict_avalanche_criterion_ok {
            return false;
        }
        let possible_reversable_bit_position = (0..Self::NUM_VARIABLES).into_iter().collect::<Vec<usize>>();
        (2..=k).into_iter().all(|criterion_degree| {
            //let possible_combinations_count = num_integer::binomial(self.input_dimension, criterion_degree);
            CombinationIterator::new(&possible_reversable_bit_position, criterion_degree)
                .all(|combination| {
                    let mut bit_mask = 0;
                    for &bit_position in combination {
                        bit_mask |= (1 << bit_position) as u32;
                    }
                    let function_equal_mask_count = (0..=max_input_value).into_iter().filter(|&x| {
                        let x_prime = x ^ bit_mask;
                        Self::compute_cellular_automata_rule(rule_number, x) == Self::compute_cellular_automata_rule(rule_number, x_prime)
                    }).count();
                    function_equal_mask_count == (1 << (Self::NUM_VARIABLES - 1))
                })
        })
    }
}