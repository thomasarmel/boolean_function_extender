use std::arch::x86_64::_popcnt32;

pub mod u512_tester;
pub mod u32_tester;

pub trait BooleanFunctionTester {
    type UnsignedRepr;

    const NUM_VARIABLES: usize;

    fn fast_bool_anf_transform_unsigned(rule_number: Self::UnsignedRepr, num_variables_function: usize) -> Self::UnsignedRepr;

    fn get_function_degree(rule_number: Self::UnsignedRepr) -> usize;

    fn is_strict_avalanche_criterion_ok(rule_number: Self::UnsignedRepr) -> bool;

    fn compute_cellular_automata_rule(rule_number: Self::UnsignedRepr, input_bits: u32) -> bool;

    fn is_function_balanced(rule_number: Self::UnsignedRepr) -> bool;

    fn is_first_order_correlation_immune(rule_number: Self::UnsignedRepr) -> bool;

    fn fast_walsh_transform(rule_number: Self::UnsignedRepr, w: u32) -> i32;

    fn is_propagation_criterion_deg_k_ok(rule_number: Self::UnsignedRepr, k: usize) -> bool;

    #[inline]
    fn fast_binary_dot_product(a: u32, b: u32) -> u32 {
        unsafe {
            _popcnt32((a & b) as i32) as u32
        }
    }
}