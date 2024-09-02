use std::collections::HashMap;
use gen_combinations::CombinationIterator;
use numext_fixed_uint::U512;
use crate::BooleanFunctionTester;

pub struct U512Tester;

impl BooleanFunctionTester for U512Tester {
    type UnsignedRepr = U512;
    const NUM_VARIABLES: usize = 9;

    const MAX_INPUT_VALUE: u32 = 2u32.pow(Self::NUM_VARIABLES as u32) - 1;
    const MAX_FUNCTION_NUMBER: Self::UnsignedRepr = U512::max_value();

    fn fast_bool_anf_transform_unsigned(rule_number: &Self::UnsignedRepr, num_variables_function: usize) -> Self::UnsignedRepr {
        const U0: U512 = U512::zero();

        let mut blocksize = 1usize;
        let mut final_f = rule_number.clone();
        for _ in 0..num_variables_function {
            let mut source = 0usize;
            while source < (1 << num_variables_function) {
                let target = source + blocksize;
                for i in 0..blocksize {
                    let f_target_i: bool = ((final_f.clone() >> (target + i)) & U512::one()) != U0;
                    let f_source_i: bool = ((final_f.clone() >> (source + i)) & U512::one()) != U0;
                    let f_target_i_xor_f_source_i = f_target_i ^ f_source_i;
                    if f_target_i_xor_f_source_i {
                        final_f = final_f | (U512::one() << (target.clone() + i));
                    } else {
                        final_f = final_f & !(U512::one() << (target.clone() + i));
                    }
                }
                source = source.clone() + (blocksize << 1);
            }
            blocksize = blocksize << 1;
        }
        final_f
    }

    fn get_function_degree(rule_number: &Self::UnsignedRepr) -> usize {
        let anf_form = Self::fast_bool_anf_transform_unsigned(rule_number, Self::NUM_VARIABLES);
        (0..=Self::MAX_INPUT_VALUE).into_iter().map(|bit_position| {
            if anf_form.clone() & (U512::one() << bit_position) != U512::zero() {
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
        (rule_number & (U512::one() << input_bits)) != U512::zero()
    }

    fn is_function_balanced(rule_number: &Self::UnsignedRepr) -> bool {
        const EXPECTED_SET_NUMBER: i32 = 1 << 8;
        rule_number.count_ones() as i32 == EXPECTED_SET_NUMBER
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
        [rule_number.clone(), Self::reverse_function(rule_number)].iter().any(|rule| {
            let mut equivalent_xor_function: U512 = U512::zero();
            for i in 0..=Self::MAX_INPUT_VALUE {
                let mut equivalent_xor_function_eval_i = false;
                for j in 0..Self::NUM_VARIABLES {
                    if rule & (U512::one() << (1 << j)) != U512::zero() {
                        equivalent_xor_function_eval_i ^= (i & (1 << j)) == 0;
                    }
                }
                equivalent_xor_function |= U512::from(equivalent_xor_function_eval_i as u32) << i;
            }
            *rule == equivalent_xor_function || *rule == Self::reverse_function(&equivalent_xor_function)
        })
    }
}

impl U512Tester {
    fn reverse_function(rule_number: &U512) -> U512 {
        !rule_number & Self::MAX_FUNCTION_NUMBER
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use numext_fixed_uint::U512;
    use crate::BooleanFunctionTester;

    #[test]
    fn test_get_function_degree() {
        assert_eq!(super::U512Tester::get_function_degree(&U512::zero()), 0);
        assert_eq!(super::U512Tester::get_function_degree(&U512::from_dec_str("13407807929942597098847186317233203207567774193617110470060917943962159749266320668916332720419656517166240687192852651607093788182615166684716733579591935").unwrap()), 4);
        assert_eq!(super::U512Tester::get_function_degree(&U512::from_dec_str("13407807929942597098847186317233203207567774193617110470605594638121991820277051203989568509643346170367187143529348261601759818256683210860310861684080895").unwrap()), 5);
        assert_eq!(super::U512Tester::get_function_degree(&U512::max_value()), 0);
    }

    #[test]
    fn test_absolute_walsh_spectrum() {
        assert_eq!(super::U512Tester::absolute_walsh_spectrum(&U512::zero()), HashMap::from([(512, 1), (0, 511)]));
        assert_eq!(super::U512Tester::absolute_walsh_spectrum(&U512::from_dec_str("13407807929942597098847186317233203207567774193617110470060917943962159749266320668916332720419656517166240687192852651607093788182615166684716733579591935").unwrap()), HashMap::from([(320, 1), (0, 490), (64, 15), (128, 6)]));
        assert_eq!(super::U512Tester::absolute_walsh_spectrum(&U512::from_dec_str("13407807929942597098847186317233203207567774193617110470605594638121991820277051203989568509643346170367187143529348261601759818256683210860310861684080895").unwrap()), HashMap::from([(352, 1), (0, 480), (32, 20), (96, 10), (160, 1)]));
        assert_eq!(super::U512Tester::absolute_walsh_spectrum(&U512::max_value()), HashMap::from([(512, 1), (0, 511)]));
    }
}