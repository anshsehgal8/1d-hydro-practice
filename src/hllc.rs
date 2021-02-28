use crate::{Conserved, Primitive};


pub struct Riemann_HLLC_t {
	pub pl: Primitive,
    pub pr: Primitive,
	pub sl: f64,
	pub sr: f64,
	pub gamma: f64,
	
}

impl Riemann_HLLC_t {
	fn fl(&self) -> Conserved { self.pl.flux(self.gamma) }
	fn fr(&self) -> Conserved { self.pr.flux(self.gamma) }
	fn ul(&self) -> Conserved { self.pl.to_cons(self.gamma)}
	fn ur(&self) -> Conserved { self.pr.to_cons(self.gamma)}
	fn d_star(&self) -> Conserved { Conserved {density: 0., momentum: 1., energy: self.s_star() } }

	fn s_star(&self) -> f64 {
		let dl = self.pl.density;
		let dr = self.pr.density;
		let vl = self.pl.velocity;
		let vr = self.pr.velocity;
		let pl = self.pl.pressure;
		let pr = self.pr.pressure;
		let sl = self.sl;
		let sr = self.sr;

		(pr - pl + dl * vl * (sl - vl) - dr * vr * (sr - vr)) / (dl * (sl - vl) - dr * (sr - vr))
	}

	fn plr(&self) -> f64 {
		let dl = self.pl.density;
		let dr = self.pr.density;
		let vl = self.pl.velocity;
		let vr = self.pr.velocity;
		let pl = self.pl.pressure;
		let pr = self.pr.pressure;

		0.5 * (pl + pr + dl * (self.sl - vl) * (self.s_star() - vl) + dr * (self.sr - vr) * (self.s_star() - vr))
	}

	fn ul_star(&self) -> Conserved {
		(self.ul() * self.sl  - self.fl() + self.d_star() * self.plr()) / (self.sl - self.s_star())
	}

	fn ur_star(&self) -> Conserved {
		(self.ur() * self.sr - self.fr() +  self.d_star() * self.plr()) / (self.sr - self.s_star())
	}


	pub fn hllc_flux(&self) -> Conserved {
		if 0. <= self.sl {
			self.fl()
		}
		else if self.sl < 0. &&  0. <= self.s_star() {
			self.fl() + (self.ul_star() - self.ul()) * self.sl 
		}
		else if self.s_star() < 0. && 0. <= self.sr {
			self.fr() +  (self.ur_star() - self.ur()) * self.sr 
		}
		else {
			self.fr()
		}
	}
}