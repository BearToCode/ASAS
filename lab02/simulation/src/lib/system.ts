export type SystemParams = {
	M: number; // mass of cart
	m: number; // mass of pendulum
	l: number; // length of pendulum
	g: number; // gravitational acceleration
};

export function nonlinearSystem(params: SystemParams) {
	const M = params.M;
	const m = params.m;
	const l = params.l;
	const g = params.g;

	const f = (x: number[], u: number, d: number) => [
		x[1], // dx
		(m * l * Math.pow(x[3], 2) * Math.sin(x[1]) + m * g * Math.sin(x[2]) * Math.cos(x[2]) + d + u) /
			(M + m * Math.pow(Math.sin(x[2]), 2)), // d2x
		x[3], // dtheta
		((((M + m) * g) / l) * Math.sin(x[2]) -
			m * Math.pow(Math.sin(x[3]), 2) * Math.cos(x[2]) * Math.sin(x[2]) +
			((d + u) * Math.cos(x[2])) / l) /
			(M + m * Math.pow(Math.sin(x[2]), 2)) // d2theta
	];

	return f;
}
