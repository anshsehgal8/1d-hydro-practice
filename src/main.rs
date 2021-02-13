use std::ops::{Add, Sub, Mul, Div};

#[derive(Copy,Clone)]
struct Conserved{
	density: f64,
	momentum_density: f64,
}

impl Add for Conserved {
	type Output = Conserved;

	fn add(self, a:Conserved) -> Conserved{
		Conserved{
			density:self.density + a.density,
			momentum_density: self.momentum_density + a.momentum_density
		}

	}

}


impl Sub for Conserved {
	type Output = Conserved;

	fn sub(self, a:Conserved) -> Conserved{
		Conserved{
			density:self.density - a.density,
			momentum_density: self.momentum_density - a.momentum_density
		}

	}

}


impl Mul<f64> for Conserved {
	type Output = Conserved;

	fn mul(self, a:f64) -> Conserved{
		Conserved{
			density:self.density * a,
			momentum_density: self.momentum_density * a
		}

	}

}


impl Div<f64> for Conserved {
	type Output = Conserved;

	fn div(self, a:f64) -> Conserved{
		Conserved{
			density:self.density / a,
			momentum_density: self.momentum_density / a
		}

	}

}


impl Conserved {
	fn velocity(self) -> f64{
		self.momentum_density / self.density
	}

	fn pressure(self, gamma:f64) -> f64{
		self.density.powf(gamma)
	}

	fn flux(self, gamma:f64) -> Conserved{
		// Conserved{
		// 	density: self.momentum_density,
		// 	momentum_density: self.density * self.velocity().powi(2) + self.pressure(gamma)
		// }

		self * self.velocity() + Conserved{ density: 0.0, momentum_density: self.pressure(gamma)}

	}

	fn sound_speed(self, gamma:f64) -> f64{
		(gamma * self.density.powf(gamma - 1.)).sqrt()
	}
}



fn hll_flux(ul: Conserved, ur: Conserved, gamma: f64) -> Conserved {
	let fl = ul.flux(gamma);
	let fr = ur.flux(gamma);
	let lambda_left_minus  = ul.velocity() - ul.sound_speed(gamma);
	let lambda_left_plus   = ul.velocity() + ul.sound_speed(gamma);
	let lambda_right_minus = ur.velocity() - ur.sound_speed(gamma);
	let lambda_right_plus  = ur.velocity() + ur.sound_speed(gamma);
	let alpha_plus  = (lambda_left_plus).max(lambda_right_plus).max(0.);
	let alpha_minus = (-lambda_left_minus).max(-lambda_right_minus).max(0.);

	((fl * alpha_plus) - (fr * alpha_minus) - (ul - ur) * alpha_plus * alpha_minus) / (alpha_plus - alpha_minus)



}





fn next(u: Vec<Conserved>, dx: f64, dt: f64, gamma: f64) -> Vec<Conserved> {
	let n = u.len();
	let mut u1 = vec![Conserved{density:0., momentum_density:0.}; n];

	for i in 1..n-1{

		let f_imh = hll_flux(u[i-1], u[i], gamma);
		let f_iph = hll_flux(u[i], u[i+1], gamma);

		u1[i] = u[i] - (f_iph - f_imh) * dt / dx
	}

	u1[0] = u[0];
	u1[n-1] = u[n-1]; 

	u1

}





fn main() {
    let x_0 = -1.;
    let x_f = 1.;
    let num_cells = 10;
    let dx = (x_f - x_0) / (num_cells as f64);
    let tfinal = 0.25;
    let dt = tfinal / 100. ; 
    let gamma = 1.;

    let xc:Vec<_> = (0..num_cells).map(|i| x_0 + (i as f64 + x_f) * dx).collect();
    let mut u: Vec<_> = xc.iter().map(|_| Conserved {density: 1., momentum_density: 0.}).collect();
    let mut t = 0.0;

    while t < tfinal {
    	u = next(u, dx, dt, gamma);
    	t += dt ;
    	println!("t = {:?}",t);

    }

}
