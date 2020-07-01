use pairing::bls12_381::{Fq, Fr, G1Affine};
use pairing::{CurveAffine, CurveProjective};
use pairing::ff::{Field, PrimeField, PrimeFieldRepr};

#[cfg(test)]
mod tests {

    use super::*;

    fn print_mont_form_value<F: PrimeField>(el: F) {
        let repr = el.into_raw_repr();
        println!("Mont repr = {}", repr);
    }

    fn add_modulus_and_print<F: PrimeField>(el: F) {
        let c = F::char();
        let mut repr = el.into_raw_repr();
        repr.add_nocarry(&c);
        println!("With extra modulus {}", repr);
    }

    fn add_modulus_multiple_and_print<F: PrimeField>(el: F, multiple: usize) {
        let c = F::char();
        let mut repr = el.into_raw_repr();
        for _ in 0..multiple {
            repr.add_nocarry(&c);
        }
        println!("With extra modulus {}", repr);
    }

    #[test]
    fn add_assign_unequal() {
        let scalar = Fr::from_str("12345").unwrap();
        let (p1_x, p1_y) = G1Affine::one().into_xy_unchecked();

        let (mut p2_x, mut p2_y) = G1Affine::one().mul(scalar.into_repr()).into_affine().into_xy_unchecked();
        let mut p2_z = Fq::one();

        // get montgomery form 

        println!("Montgomery form of P2.X");
        print_mont_form_value(p2_x);

        // we add p1 (affine point) to p2 (projective)

        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl

        // Z1Z1 = Z1^2
        let mut z1z1 = p2_z;
        z1z1.square();

        // U2 = X2*Z1Z1
        let mut u2 = p1_x;
        u2.mul_assign(&z1z1);

        // S2 = Y2*Z1*Z1Z1
        let mut s2 = p1_y;
        s2.mul_assign(&p2_z);
        s2.mul_assign(&z1z1);

        if p2_x == u2 && p2_y == s2 {
            panic!("we add unequal points");
            // The two points are equal, so we double.
        } else {
            // If we're adding -a and a together, p2_z becomes zero as H becomes zero.

            // H = U2-X1
            let mut h = u2;
            h.sub_assign(&p2_x);

            // HH = H^2
            let mut hh = h;
            hh.square();

            // I = 4*HH
            let mut i = hh;
            i.double();
            i.double();

            // here is how we do to get partial reduction
            let mut i_repr = hh.into_raw_repr(); // into_raw_repr returns montgomery form (but fully reduced!);
            i_repr.shl(2); // shift left by 2 bits

            println!("I unreduced value");
            println!("4*HH = {}", i_repr);

            println!("I reduced value");
            print_mont_form_value(i);

            println!("See equality");
            for multiple in 0..=4 {
                println!("+ {} moduluses", multiple);
                add_modulus_multiple_and_print(i, multiple);
            }

            // J = H*I
            let mut j = h;
            j.mul_assign(&i);

            // r = 2*(S2-Y1)
            let mut r = s2;
            r.sub_assign(&p2_y);
            r.double();

            // V = X1*I
            let mut v = p2_x;
            v.mul_assign(&i);

            println!("Fully reduced V in montgomery form = {}", v.into_raw_repr());

            // V can be not fully reduced and have a value 
            add_modulus_and_print(v);

            // X3 = r^2 - J - 2*V
            p2_x = r;
            p2_x.square();
            p2_x.sub_assign(&j);
            p2_x.sub_assign(&v);
            p2_x.sub_assign(&v);

            // Y3 = r*(V-X3)-2*Y1*J
            j.mul_assign(&p2_y); // J = 2*Y1*J
            j.double();
            p2_y = v;
            p2_y.sub_assign(&p2_x);
            p2_y.mul_assign(&r);
            p2_y.sub_assign(&j);

            // Z3 = (Z1+H)^2-Z1Z1-HH
            p2_z.add_assign(&h);
            p2_z.square();
            p2_z.sub_assign(&z1z1);
            p2_z.sub_assign(&hh);
        }

        println!("Result in normal form");
        println!("X = {}, Y = {}, Z = {}", p2_x, p2_y, p2_z);
    }
}


