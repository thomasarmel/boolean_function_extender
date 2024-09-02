use std::collections::HashMap;
use gen_combinations::CombinationIterator;
use crate::BooleanFunctionTester;

pub struct U32Tester;

impl BooleanFunctionTester for U32Tester {
    type UnsignedRepr = u32;
    const NUM_VARIABLES: usize = 5;
    const MAX_INPUT_VALUE: u32 = 2u32.pow(Self::NUM_VARIABLES as u32) - 1;
    const MAX_FUNCTION_NUMBER: Self::UnsignedRepr = u32::MAX;

    fn fast_bool_anf_transform_unsigned(rule_number: &Self::UnsignedRepr, num_variables_function: usize) -> Self::UnsignedRepr {
        fast_boolean_anf_transform::fast_bool_anf_transform_unsigned(*rule_number, num_variables_function)
    }

    fn get_function_degree(rule_number: &Self::UnsignedRepr) -> usize {
        let anf_form = Self::fast_bool_anf_transform_unsigned(rule_number, Self::NUM_VARIABLES);
        (0..=Self::MAX_INPUT_VALUE).into_iter().map(|bit_position| {
            if anf_form & (1 << bit_position) != 0 {
                bit_position.count_ones() as usize
            } else {
                0
            }
        }).max().unwrap_or(0)
    }

    fn is_strict_avalanche_criterion_ok(rule_number: &Self::UnsignedRepr) -> bool {
        (0..Self::NUM_VARIABLES).into_iter().all(|constant_position| {
            let constant = 1 << constant_position;
            (0..=Self::MAX_INPUT_VALUE).into_iter().filter(|&x| {
                let x_prime = x ^ constant;
                Self::compute_cellular_automata_rule(rule_number, x) == Self::compute_cellular_automata_rule(rule_number, x_prime)
            }).count() == (1 << (Self::NUM_VARIABLES - 1))
        })
    }

    fn compute_cellular_automata_rule(rule_number: &Self::UnsignedRepr, input_bits: u32) -> bool {
        #[cfg(debug_assertions)]
        if input_bits > Self::MAX_INPUT_VALUE {
            panic!("Input bits must be less or equal than {}", Self::MAX_INPUT_VALUE);
        }
        (rule_number & (1 << input_bits)) != 0
    }

    fn is_function_balanced(rule_number: &Self::UnsignedRepr) -> bool {
        const EXPECTED_SET_NUMBER: u32 = 1 << 4;
        rule_number.count_ones() == EXPECTED_SET_NUMBER
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

    // https://doc.sagemath.org/html/en/reference/cryptography/sage/crypto/boolean_function.html#sage.crypto.boolean_function.BooleanFunction.walsh_hadamard_transform
    fn fast_walsh_transform(rule_number: &Self::UnsignedRepr, w: u32) -> i32 {
        (0..=Self::MAX_INPUT_VALUE).map(|x| {
            if (Self::compute_cellular_automata_rule(rule_number, x) as u32 + Self::fast_binary_dot_product(w, x as u32)) & 1 == 0 { // % modulo 2
                1
            } else {
                -1
            }
        }).sum()
    }

    fn absolute_walsh_spectrum(rule_number: &Self::UnsignedRepr) -> HashMap<u32, usize> {
        let mut absolute_walsh_value_count_map: HashMap<u32, usize> = HashMap::new();
        (0..=Self::MAX_INPUT_VALUE)
            .for_each(|w| {
                let absolute_walsh_value = Self::fast_walsh_transform(rule_number, w).unsigned_abs();
                if !absolute_walsh_value_count_map.contains_key(&absolute_walsh_value) {
                    absolute_walsh_value_count_map.insert(absolute_walsh_value, 1);
                } else {
                    let count = absolute_walsh_value_count_map.get_mut(&absolute_walsh_value).unwrap();
                    *count += 1;
                }
            });
        absolute_walsh_value_count_map
    }

    fn is_propagation_criterion_deg_k_ok(rule_number: &Self::UnsignedRepr, k: usize) -> bool {
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
                    let function_equal_mask_count = (0..=Self::MAX_INPUT_VALUE).into_iter().filter(|&x| {
                        let x_prime = x ^ bit_mask;
                        Self::compute_cellular_automata_rule(rule_number, x) == Self::compute_cellular_automata_rule(rule_number, x_prime)
                    }).count();
                    function_equal_mask_count == (1 << (Self::NUM_VARIABLES - 1))
                })
        })
    }

    fn is_function_linear(rule_number: &Self::UnsignedRepr) -> bool {
        [*rule_number, Self::reverse_function(*rule_number)].iter().any(|rule| {
            let mut equivalent_xor_function: u32 = 0;
            for i in 0..=Self::MAX_INPUT_VALUE {
                let mut equivalent_xor_function_eval_i = false;
                for j in 0..Self::NUM_VARIABLES {
                    if *rule & (1 << (1 << j)) != 0 {
                        equivalent_xor_function_eval_i ^= (i & (1 << j)) == 0;
                    }
                }
                equivalent_xor_function |= (equivalent_xor_function_eval_i as u32) << i;
            }
            *rule == equivalent_xor_function || *rule == Self::reverse_function(equivalent_xor_function)
        })
    }
}

impl U32Tester {
    pub fn fast_auto_correlation_transform(rule_number: &u32, w: u32) -> i32 {
        (0..=Self::MAX_INPUT_VALUE).map(|x| {
            if Self::compute_cellular_automata_rule(rule_number, x) ^ Self::compute_cellular_automata_rule(rule_number, x ^ w) {
                -1
            } else {
                1
            }
        }).sum()
    }

    pub fn absolute_autocorrelation_spectrum(rule_number: &u32) -> HashMap<u32, usize> {
        let mut absolute_autocorrelation_value_count_map: HashMap<u32, usize> = HashMap::new();
        (0..=Self::MAX_INPUT_VALUE)
            .for_each(|w| {
                let absolute_autocorrelation_value = Self::fast_auto_correlation_transform(rule_number, w).unsigned_abs();
                if !absolute_autocorrelation_value_count_map.contains_key(&absolute_autocorrelation_value) {
                    absolute_autocorrelation_value_count_map.insert(absolute_autocorrelation_value, 1);
                } else {
                    let count = absolute_autocorrelation_value_count_map.get_mut(&absolute_autocorrelation_value).unwrap();
                    *count += 1;
                }
            });
        absolute_autocorrelation_value_count_map
    }

    fn reverse_function(rule_number: u32) -> u32 {
        !rule_number & Self::MAX_FUNCTION_NUMBER
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use crate::BooleanFunctionTester;

    #[test]
    fn test_get_function_degree() {
        assert_eq!(super::U32Tester::get_function_degree(&0), 0);
        assert_eq!(super::U32Tester::get_function_degree(&3755921403), 4);
        assert_eq!(super::U32Tester::get_function_degree(&3755921407), 5);
        assert_eq!(super::U32Tester::get_function_degree(&4294967295), 0);
    }

    #[test]
    fn test_absolute_walsh_spectrum() {
        assert_eq!(super::U32Tester::absolute_walsh_spectrum(&0), HashMap::from([(32, 1), (0, 31)]));
        assert_eq!(super::U32Tester::absolute_walsh_spectrum(&3755921403), HashMap::from([(20, 1), (0, 10), (8, 6), (4, 15)]));
        assert_eq!(super::U32Tester::absolute_walsh_spectrum(&3755921407), HashMap::from([(22, 1), (2, 20), (10, 1), (6, 10)]));
        assert_eq!(super::U32Tester::absolute_walsh_spectrum(&4294967295), HashMap::from([(32, 1), (0, 31)]));
    }

    #[test]
    fn test_absolute_autocorrelation_spectrum() {
        assert_eq!(super::U32Tester::absolute_autocorrelation_spectrum(&0), HashMap::from([(32, 32)]));
        assert_eq!(super::U32Tester::absolute_autocorrelation_spectrum(&3755921403), HashMap::from([(32, 1), (16, 15), (8, 16)]));
        assert_eq!(super::U32Tester::absolute_autocorrelation_spectrum(&3755921407), HashMap::from([(32, 1), (20, 10), (12, 21)]));
        assert_eq!(super::U32Tester::absolute_autocorrelation_spectrum(&4294967295), HashMap::from([(32, 32)]));
    }
}