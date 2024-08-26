use std::arch::x86_64::_bzhi_u32;
use gen_combinations::CombinationIterator;
use numext_fixed_uint::U512;
use crate::BooleanFunctionTester;

pub struct U512Tester;

impl BooleanFunctionTester for U512Tester {
    type UnsignedRepr = U512;
    const NUM_VARIABLES: usize = 9;

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
        let max_input_value = unsafe { _bzhi_u32(u32::MAX, Self::NUM_VARIABLES as u32) };
        (0..=max_input_value).into_iter().map(|bit_position| {
            if anf_form.clone() & (U512::one() << bit_position) != U512::zero() {
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
        #[cfg(debug_assertions)]
        {
            let max_input_value = unsafe { _bzhi_u32(u32::MAX, Self::NUM_VARIABLES as u32) };
            if input_bits > max_input_value {
                panic!("Input bits must be less or equal than {}", max_input_value);
            }
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