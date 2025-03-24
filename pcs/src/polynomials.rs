use blstrs::Scalar;
use num_traits::pow;
use rand::RngCore;

use core::ops::Div;
use group::ff::Field;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Polynomial(pub Vec<Scalar>);

impl Polynomial {
    pub fn new(scalars: &[Scalar]) -> Self {
        Polynomial(scalars.to_vec())
    }

    pub fn new_from_bytes(bytes: &[u8]) -> Self {
        let scalars: Vec<Scalar> = bytes.into_iter().map(|d| Scalar::from(d.clone() as u64)).collect();
        Polynomial(scalars.to_vec())
    }

    pub fn evaluate(&self, point: Scalar) -> Scalar {
        let mut total = Scalar::from(0_u64);
        let mut cur_p = Scalar::from(1_u64);
        for (_, cof) in self.0.iter().enumerate() {
            total += cur_p * cof;
            cur_p = cur_p * point;
        }
        total
    }

    fn is_zero(&self) -> bool {
        self.0.is_empty() || self.0.iter().all(|a| a.is_zero().into())
    }

    fn leading_coefficient(&self) -> Option<Scalar> {
        self.0.last().copied()
    }
}

impl Div for Polynomial {
    type Output = Self;
    fn div(self, divisor: Self) -> Self::Output {
        if self.is_zero() {
            panic!("Dividing by zero poly!")
        } else if self.0.len() < divisor.0.len() {
            Polynomial::new(&[Scalar::from(0)])
        } else {
            let mut quotient = Polynomial::new(&vec![Scalar::ZERO; self.0.len() - divisor.0.len() + 1]);
            let mut remainder: Polynomial = self.clone().into();

            let divisor_leading_inv = divisor.leading_coefficient().unwrap().invert().unwrap();
            while !remainder.is_zero() && remainder.0.len() >= divisor.0.len() {
                let cur_q_coeff = remainder.leading_coefficient().unwrap() * divisor_leading_inv;
                let cur_q_degree = remainder.0.len() - divisor.0.len();
                quotient.0[cur_q_degree] = cur_q_coeff;

                for (i, div_coeff) in divisor.0.iter().enumerate() {
                    remainder.0[cur_q_degree + i] -= &(cur_q_coeff * div_coeff);
                }
                while let Some(true) = remainder.0.last().map(|c| c.is_zero().into()) {
                    remainder.0.pop();
                }
            }
            
            quotient
        }
    }
}