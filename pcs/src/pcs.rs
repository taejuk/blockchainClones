use crate::polynomials::Polynomial;
use blstrs::{pairing, G1Projective, G2Projective, Scalar};
use group::Curve;
use group::{ff::Field as FieldT, Group};
use rand::Rng;

fn generate_tau_points<T: Group + std::ops::Mul<Scalar, Output=T>>(
    generator: T,
    tau: Scalar,
    length: usize
) -> Vec<T> {
    let mut generators = Vec::with_capacity(length);
    generators.push(generator);
    let mut g = generator.clone();

    for _ in 1..length {
        g = g * tau;
        generators.push(g);
    }

    generators
}

#[derive(Clone, Debug)]
pub struct GlobalParameters {
    pub gs: Vec<G1Projective>,
    hs: Vec<G2Projective>
}

impl GlobalParameters {
    fn new(gs: Vec<G1Projective>, hs: Vec<G2Projective>) -> Self {
        GlobalParameters {gs, hs}
    }
}

#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    IncorrectDegree,
    SetupIncomplete,
}

pub trait PolynomialCommitment {
    fn setup(&mut self, d: usize) -> GlobalParameters;

    fn commit(&self, polynomial: &Polynomial) -> Result<G1Projective, Error>;
    fn create_witness(&self, polynomial: Polynomial, point: Scalar) -> (G1Projective, Scalar);
    fn verify_evaluation(
        &self,
        committed_polynomial: G1Projective,
        point: Scalar,
        evaluation: Scalar,
        witness: G1Projective,
    ) -> bool;
}

#[derive(Debug)]
pub struct GenericPolynomialCommitment {
    global_parameters: Option<GlobalParameters>,
}

impl GenericPolynomialCommitment {
    // This might seem useless for now. I am keeping it, as I might want to come back later for more initialization values
    pub fn new() -> Self {
        GenericPolynomialCommitment {
            global_parameters: None,
        }
    }
}

impl PolynomialCommitment for GenericPolynomialCommitment {
    fn setup(
        &mut self,
        d: usize
    ) -> GlobalParameters {
        let mut rng = rand::rng();
        let tau: u64 = rng.random();
        let tau = Scalar::from(tau);
        let gs = generate_tau_points(G1Projective::generator(), tau, d);
        let hs = generate_tau_points(G2Projective::generator(), tau, d);
        let globalparameters = GlobalParameters::new(gs, hs);
        self.global_parameters = Some(globalparameters.clone());
        globalparameters
    }

    fn commit(&self, polymial: &Polynomial) -> Result<G1Projective, Error> {
        if self.global_parameters.is_none() {
            return Err(Error::SetupIncomplete);
        }

        let global_parameters = &self.global_parameters.as_ref().unwrap();
        if polymial.0.len() != global_parameters.gs.len() {
            return Err(Error::IncorrectDegree);
        }

        Ok(G1Projective::multi_exp(&global_parameters.gs, &polymial.0))
    }

    fn create_witness(&self, polynomial: Polynomial, point: Scalar) -> (G1Projective, Scalar) {
        let evaluation = polynomial.evaluate(point);
        let mut witness_polynomial = polynomial.clone();
        witness_polynomial.0[0] -= &evaluation;
        let divisor = Polynomial::new(&[-point, Scalar::ONE]);
        witness_polynomial = witness_polynomial / divisor;
        let witness = G1Projective::multi_exp(
            &self.global_parameters.as_ref().unwrap().gs[..witness_polynomial.0.len()],
            &witness_polynomial.0
        );

        (witness, evaluation)
    }

    fn verify_evaluation(
        &self,
        committed_polynomial: G1Projective,
        point: Scalar,
        evaluation: Scalar,
        witness: G1Projective,
    ) -> bool {
        let g1 = G1Projective::generator();
        let g2 = G2Projective::generator();
        let evaluation_inverse = g1 * -evaluation;

        let left_pairing = committed_polynomial + evaluation_inverse;
        let lhs = pairing(&left_pairing.to_affine(), &g2.to_affine());

        let point_commitment_inverted = g2 * -point;

        let right_side = self.global_parameters.as_ref().unwrap().hs[1] + point_commitment_inverted;
        let rhs = pairing(&witness.to_affine(), &right_side.to_affine());
        lhs == rhs
    }

}