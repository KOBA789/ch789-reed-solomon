use std::{
    fmt::{self, Debug, Display},
    ops::{Add, Mul},
};

use once_cell::sync::Lazy;

#[derive(Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Hash)]
struct AlphaPoly(u8);

impl Debug for AlphaPoly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        //write!(f, "{:#010b}", self.0)
        if self.0 == 0 {
            write!(f, "0")
        } else {
            write!(f, "a^{}", ALPHA_REV_TABLE[self.0 as usize])
        }
    }
}

static ALPHA_TABLE: Lazy<[AlphaPoly; 255]> = Lazy::new(AlphaPoly::generate_table);

static ALPHA_REV_TABLE: Lazy<[u8; 256]> = Lazy::new(|| {
    let mut rev_table = [0u8; 256];
    for (idx, alpha) in ALPHA_TABLE.iter().enumerate() {
        rev_table[alpha.0 as usize] = idx as u8;
    }
    rev_table
});

impl AlphaPoly {
    const GENERATOR: AlphaPoly = AlphaPoly(0b11101);

    fn mul_alpha(&self) -> AlphaPoly {
        if self.0 & 0b1000_0000 == 0 {
            AlphaPoly(self.0 << 1)
        } else {
            AlphaPoly(self.0 << 1) + Self::GENERATOR
        }
    }

    fn generate_table() -> [AlphaPoly; 255] {
        let mut alpha_table = [AlphaPoly(0); 255];
        alpha_table[0] = AlphaPoly(1);
        for i in 1..u8::MAX {
            alpha_table[i as usize] = alpha_table[i as usize - 1].mul_alpha();
        }
        alpha_table
    }

    fn inv(self) -> AlphaPoly {
        let idx = ALPHA_REV_TABLE[self.0 as usize];
        if idx == 0 {
            return ALPHA_TABLE[0];
        }
        let idx = 255 - idx;
        ALPHA_TABLE[idx as usize]
    }
}

impl Add for AlphaPoly {
    type Output = AlphaPoly;

    fn add(self, rhs: Self) -> Self::Output {
        #[allow(clippy::suspicious_arithmetic_impl)]
        AlphaPoly(self.0 ^ rhs.0)
    }
}

impl Mul for AlphaPoly {
    type Output = AlphaPoly;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.0 == 0 || rhs.0 == 0 {
            return AlphaPoly(0);
        }
        let lhs_idx = ALPHA_REV_TABLE[self.0 as usize];
        let rhs_idx = ALPHA_REV_TABLE[rhs.0 as usize];
        let mul_idx = (lhs_idx as u16 + rhs_idx as u16) % 255;
        ALPHA_TABLE[mul_idx as usize]
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Poly(Vec<AlphaPoly>);

impl Display for Poly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let bytes = self
            .0
            .iter()
            .map(|alpha| ALPHA_REV_TABLE[alpha.0 as usize])
            .collect::<Vec<u8>>();
        write!(f, "{}", bytes.escape_ascii())
    }
}

impl Poly {
    fn normalize(&mut self) {
        let trailing_zeros = self.0.iter().rev().take_while(|alpha| alpha.0 == 0).count();
        self.0.truncate(self.0.len() - trailing_zeros);
    }
}

impl Add for Poly {
    type Output = Poly;

    fn add(mut self, mut rhs: Self) -> Self::Output {
        let new_len = self.0.len().max(rhs.0.len());
        self.0.resize(new_len, AlphaPoly(0));
        rhs.0.resize(new_len, AlphaPoly(0));
        for (lhs, rhs) in self.0.iter_mut().zip(rhs.0.into_iter()) {
            *lhs = *lhs + rhs;
        }
        self.normalize();
        self
    }
}

impl Mul for Poly {
    type Output = Poly;

    fn mul(self, rhs: Self) -> Self::Output {
        let len = self.0.len() + rhs.0.len() - 1;
        let mut ret = vec![AlphaPoly(0); len];
        for (lhs_x, lhs_alpha) in self.0.iter().enumerate() {
            for (rhs_x, rhs_alpha) in rhs.0.iter().enumerate() {
                ret[lhs_x + rhs_x] = ret[lhs_x + rhs_x] + (*lhs_alpha * *rhs_alpha);
            }
        }
        let mut poly = Poly(ret);
        poly.normalize();
        poly
    }
}

impl Mul<AlphaPoly> for Poly {
    type Output = Poly;

    fn mul(self, rhs: AlphaPoly) -> Self::Output {
        self * Poly(vec![rhs])
    }
}

impl Poly {
    fn div_mod(self, rhs: Self) -> (Poly, Poly) {
        let mut lhs = self;
        let mut q = Poly(vec![]);
        loop {
            if lhs.0.len() < rhs.0.len() {
                return (q, lhs);
            }
            let x_idx = lhs.0.len() - rhs.0.len();
            let mut p = Poly(vec![AlphaPoly(0); x_idx + 1]);
            let lhs_alpha_idx = ALPHA_REV_TABLE[lhs.0[lhs.0.len() - 1].0 as usize] as i16;
            let rhs_alpha_idx = ALPHA_REV_TABLE[rhs.0[rhs.0.len() - 1].0 as usize] as i16;
            let mut alpha_idx = lhs_alpha_idx - rhs_alpha_idx;
            if alpha_idx < 0 {
                alpha_idx += 255;
            }
            let len = p.0.len();
            p.0[len - 1] = ALPHA_TABLE[alpha_idx as usize];
            let m = rhs.clone() * p.clone() + lhs;
            lhs = m;
            q = q + p;
        }
    }

    fn assign_alpha(self, alpha: AlphaPoly) -> AlphaPoly {
        let (cum, _alpha_exp) = self
            .0
            .iter()
            .fold((AlphaPoly(0), ALPHA_TABLE[0]), |(cum, alpha_exp), coeff| {
                (cum + *coeff * alpha_exp, alpha_exp * alpha)
            });
        cum
    }

    fn gcd(self, rhs: Self, t: usize, a_2: Poly, a_1: Poly) -> (Poly, Poly) {
        let lhs = self;
        if rhs.0.len() <= t {
            return (rhs, a_1);
        }
        let (q, rem) = lhs.div_mod(rhs.clone());
        let a = a_2 + q * a_1.clone();
        rhs.gcd(rem, t, a_1, a)
    }
}

fn main() {}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    #[test]
    fn test_alpha_table_uniqueness() {
        let table = AlphaPoly::generate_table();
        let count = table.iter().sorted().unique().count();
        assert_eq!(count, 255);
    }

    #[test]
    fn test_alpha_poly_mul() {
        assert_eq!(ALPHA_TABLE[6] * ALPHA_TABLE[7], ALPHA_TABLE[13]);
        assert_eq!(ALPHA_TABLE[200] * ALPHA_TABLE[100], ALPHA_TABLE[45]);
        assert_eq!(ALPHA_TABLE[5] * AlphaPoly(0), AlphaPoly(0));
    }

    #[test]
    fn test_poly_add() {
        assert_eq!(
            Poly(vec![AlphaPoly(0)]) + Poly(vec![AlphaPoly(0)]),
            Poly(vec![])
        );
        assert_eq!(
            Poly(vec![ALPHA_TABLE[0]]) + Poly(vec![AlphaPoly(0)]),
            Poly(vec![ALPHA_TABLE[0]])
        );
        assert_eq!(
            Poly(vec![AlphaPoly(0), ALPHA_TABLE[1]]) + Poly(vec![ALPHA_TABLE[2]]),
            Poly(vec![ALPHA_TABLE[2], ALPHA_TABLE[1]])
        );
    }

    #[test]
    fn test_normalize() {
        let mut poly = Poly(vec![
            AlphaPoly(0),
            ALPHA_TABLE[1],
            AlphaPoly(0),
            AlphaPoly(0),
        ]);
        poly.normalize();
        assert_eq!(poly, Poly(vec![AlphaPoly(0), ALPHA_TABLE[1]]));
    }

    #[test]
    fn test_poly_mul() {
        assert_eq!(
            Poly(vec![
                AlphaPoly(0),
                AlphaPoly(0),
                ALPHA_TABLE[1],
                ALPHA_TABLE[2]
            ]) * Poly(vec![ALPHA_TABLE[1], ALPHA_TABLE[1]]),
            Poly(vec![
                AlphaPoly(0),
                AlphaPoly(0),
                ALPHA_TABLE[2],
                ALPHA_TABLE[3] + ALPHA_TABLE[2],
                ALPHA_TABLE[3]
            ])
        );
    }

    #[test]
    fn test_div_mod() {
        let lhs = Poly(vec![ALPHA_TABLE[1], ALPHA_TABLE[8], ALPHA_TABLE[100]]);
        let rhs = Poly(vec![ALPHA_TABLE[0], ALPHA_TABLE[5]]);
        assert_eq!(
            lhs.div_mod(rhs),
            (
                Poly(vec![ALPHA_TABLE[170], ALPHA_TABLE[95]]),
                Poly(vec![ALPHA_TABLE[157]])
            ),
        );

        let lhs = Poly(vec![ALPHA_TABLE[1], ALPHA_TABLE[8], ALPHA_TABLE[5]]);
        let rhs = Poly(vec![ALPHA_TABLE[0], ALPHA_TABLE[100]]);
        assert_eq!(
            lhs.div_mod(rhs),
            (
                Poly(vec![ALPHA_TABLE[134], ALPHA_TABLE[160]]),
                Poly(vec![ALPHA_TABLE[251]])
            ),
        );

        let lhs = Poly(vec![ALPHA_TABLE[1], ALPHA_TABLE[8], ALPHA_TABLE[5]]);
        let rhs = Poly(vec![ALPHA_TABLE[0], AlphaPoly(0), ALPHA_TABLE[100]]);
        assert_eq!(
            lhs.div_mod(rhs),
            (
                Poly(vec![ALPHA_TABLE[160]]),
                Poly(vec![ALPHA_TABLE[156], ALPHA_TABLE[8]])
            ),
        );
    }

    #[test]
    fn test_poly_assign_alpha() {
        assert_eq!(
            Poly(vec![ALPHA_TABLE[2], ALPHA_TABLE[8], ALPHA_TABLE[100]])
                .assign_alpha(ALPHA_TABLE[2]),
            ALPHA_TABLE[243]
        );
        assert_eq!(
            Poly(vec![ALPHA_TABLE[2], ALPHA_TABLE[8], ALPHA_TABLE[100]])
                .assign_alpha(ALPHA_TABLE[200]),
            ALPHA_TABLE[71]
        );
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_reedsolomon() {
        let I = Poly(b"hello world".iter().map(|b| ALPHA_TABLE[*b as usize]).collect());
        let k = I.0.len();
        let t = 6;
        let n = k + 2 * t;
        let mut x_n_k = vec![AlphaPoly(0); n - k + 1];
        x_n_k[n - k] = ALPHA_TABLE[0];
        let x_n_k = Poly(x_n_k);
        let mut G = Poly(vec![ALPHA_TABLE[0]]);
        for idx in 1..=(2 * t) {
            G = G * Poly(vec![ALPHA_TABLE[idx], ALPHA_TABLE[0]]);
        }
        let x_n_k_I = x_n_k * I;
        let (_, P) = x_n_k_I.clone().div_mod(G);
        let C = x_n_k_I + P;
        let E = Poly(vec![
            AlphaPoly(0),
            AlphaPoly(0),
            AlphaPoly(0),
            ALPHA_TABLE[102],
            AlphaPoly(0),
            ALPHA_TABLE[102],
            AlphaPoly(0),
            AlphaPoly(0),
            ALPHA_TABLE[0],
            AlphaPoly(0),
            AlphaPoly(0),
            ALPHA_TABLE[200],
            ALPHA_TABLE[60],
            AlphaPoly(0),
            ALPHA_TABLE[60],
            AlphaPoly(0),
            AlphaPoly(0),
        ]);
        let y = C + E;
        println!("y = {}", &y);
        let mut S = Poly(vec![]);
        for idx in 1..=(2 * t) {
            let s = y.clone().assign_alpha(ALPHA_TABLE[idx]);
            S.0.push(s);
        }
        S.normalize();
        let mut x_2t = Poly(vec![AlphaPoly(0); 2 * t + 1]);
        x_2t.0[2 * t] = ALPHA_TABLE[0];
        if S == Poly(vec![]) {
            println!("no errors");
            return;
        }
        let (r, a) = x_2t.gcd(S, t, Poly(vec![]), Poly(vec![ALPHA_TABLE[0]]));
        let a_inv = a.0[0].inv();
        let eta = r * a_inv;
        let sigma = a * a_inv;
        let mut eps = vec![];
        for i in 0..n {
            let x = if i == 0 {
                ALPHA_TABLE[0]
            } else {
                ALPHA_TABLE[i].inv()
            };
            let s = sigma.clone().assign_alpha(x);
            if s == AlphaPoly(0) {
                eps.push(i);
            }
            //eprintln!("{i} {:?}", s);
        }
        let sigma2 = eps.iter().fold(Poly(vec![AlphaPoly(0)]), |sigma2, l| {
            let t = eps
                .iter()
                .filter(|i| *i != l)
                .fold(Poly(vec![ALPHA_TABLE[*l]]), |t, i| {
                    t * Poly(vec![ALPHA_TABLE[0], ALPHA_TABLE[*i]])
                });
            sigma2 + t
        });
        //dbg!(&sigma2);
        let mut e = Poly(vec![AlphaPoly(0); n]);
        for i in eps.iter() {
            let alpha_inv = ALPHA_TABLE[*i].inv();
            let numera = eta.clone().assign_alpha(alpha_inv);
            let denomi = sigma2.clone().assign_alpha(alpha_inv);
            e.0[*i] = numera * denomi.inv();
        }
        //dbg!(&e);
        let c = y + e;
        println!("c = {}", c);
    }

    #[test]
    #[ignore]
    fn test_gcd() {
        let rhs = Poly(vec![ALPHA_TABLE[0]]);
        let lhs = Poly(vec![ALPHA_TABLE[0]]);
        //assert_eq!(rhs.gcd(lhs), Poly(vec![ALPHA_TABLE[0]]));

        let rhs = Poly(vec![ALPHA_TABLE[100], ALPHA_TABLE[3]])
            * Poly(vec![ALPHA_TABLE[2], ALPHA_TABLE[0]]);
        let lhs = Poly(vec![ALPHA_TABLE[100], ALPHA_TABLE[3]])
            * Poly(vec![ALPHA_TABLE[3], ALPHA_TABLE[0]]);
        //assert_eq!(rhs.gcd(lhs), Poly(vec![ALPHA_TABLE[100], ALPHA_TABLE[3]]));
    }
}
