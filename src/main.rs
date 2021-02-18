use std::fs;
use std::io::Write;
use derive_more::{Add,Sub, Mul, Div};




/**
 * Conserved quantities
 */
#[derive(Copy,Clone, Add, Sub, Mul, Div)]
struct Conserved{
    density: f64,
    momentum: f64,
    energy: f64,
}


#[derive(Copy,Clone, Add, Sub, Mul, Div)]
struct Primitive{
    density: f64,
    velocity: f64,
    pressure: f64,
}





// ============================================================================
impl Conserved {

	fn to_prim(self, gamma:f64) -> Primitive {

		Primitive{
			density: self.density,
			velocity: self.momentum / self.density,
			pressure: (gamma - 1.) * (self. energy - 0.5 * self.momentum.powi(2) / self.density)		
		}
	}


    // fn flux(self, gamma:f64) -> Conserved {
    //     self * self.velocity() + Conserved { density: 0.0, momentum: self.pressure(gamma) }
    // }

    // fn sound_speed(self, gamma:f64) -> f64{
    //     (gamma * self.density.powf(gamma - 1.0)).sqrt()
    // }
}



impl Primitive {

	fn to_cons(self, gamma:f64) -> Conserved {

		Conserved{
			density: self.density,
			momentum: self.density * self.velocity,
			energy: self.pressure / (gamma - 1.) + 0.5 * self.density * self.velocity.powi(2)
		}

	}

	fn flux(self, gamma:f64) -> Conserved {

		self.to_cons(gamma) * self.velocity + Conserved{ density: 0., momentum: self.pressure, energy: self.pressure * self.velocity}

	}

	fn sound_speed(self, gamma:f64) -> f64{
        (gamma * self.pressure / self.density).sqrt()
    }
}



// ============================================================================
fn hll_flux(ul: Conserved, ur: Conserved, gamma: f64) -> Conserved {
	let pl = ul.to_prim(gamma);
	let pr = ur.to_prim(gamma);
    let fl = pl.flux(gamma);
    let fr = pr.flux(gamma);
    let lambda_left_minus  = pl.velocity - pl.sound_speed(gamma);
    let lambda_left_plus   = pl.velocity + pl.sound_speed(gamma);
    let lambda_right_minus = pr.velocity - pr.sound_speed(gamma);
    let lambda_right_plus  = pr.velocity + pr.sound_speed(gamma);
    let alpha_plus  = ( lambda_left_plus) .max( lambda_right_plus). max(0.0);
    let alpha_minus = (-lambda_left_minus).max(-lambda_right_minus).max(0.0);
    ((fl * alpha_plus) + (fr * alpha_minus) - (ur - ul) * alpha_plus * alpha_minus) / (alpha_plus + alpha_minus)
}




// ============================================================================
fn next(u: Vec<Conserved>, dx: f64, dt: f64, gamma: f64) -> Vec<Conserved> {

    let n = u.len();
    let mut u1 = vec![Conserved{density:0., momentum:0., energy:0.}; n];

    for i in 1..n-1 {
        let f_imh = hll_flux(u[i-1], u[i], gamma);
        let f_iph = hll_flux(u[i], u[i+1], gamma);
        u1[i] = u[i] - (f_iph - f_imh) * dt / dx
    }

    u1[0] = u[0];
    u1[n-1] = u[n-1]; 
    u1
}




// ============================================================================
fn shocktube(x: f64, x_split: f64, rho_p_left: (f64,f64), rho_p_right: (f64,f64)) -> Primitive {
    if x < x_split {
        Primitive {
            density: rho_p_left.0,
            velocity: 0.0,
            pressure: rho_p_left.1,
        }
    }
    else {
        Primitive {
            density: rho_p_right.0,
            velocity: 0.0,
            pressure: rho_p_right.1,
        }
    }
}




// ============================================================================
fn main() {

    let num_cells = 1000;
    let tfinal = 0.25;
    let x_0 = -1.0;
    let x_f =  1.0;
    let dx = (x_f - x_0) / (num_cells as f64);
    let dt = tfinal / 1000.0; 
    let gamma = 5. / 3.;

    let xc: Vec<_> = (0..num_cells).map(|i| x_0 + (i as f64 + 0.5) * dx).collect();
    let mut u: Vec<_> = xc.iter().map(|&x| shocktube(x, 0.0, (1.0, 1.0), (0.1, 0.125)).to_cons(gamma)).collect();
    let mut t = 0.0;

    while t < tfinal {
        u = next(u, dx, dt, gamma);
        t += dt ;
    }

    let file = fs::File::create("solution.dat").unwrap();

    for (xc, u) in xc.iter().zip(u) {
    	let p = u.to_prim(gamma);
        writeln!(&file, "{:.6} {:.6} {:.6} {:.6}", xc, p.density, p.velocity, p.pressure).unwrap();
    }
}
