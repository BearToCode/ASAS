import { rk4 } from './rk4';
import { scalarMul, sumArrays } from './utils';

export function ab4(f: (t: number, x: number[]) => number[], x0: number[]) {
	let t = 0;
	let prev_x = x0;

	let previousXDot: number[][] = new Array(3).fill(f(t, x0)); // Store the last 3 derivatives

	// To compute the initial steps
	const initialSolver = rk4(f, x0);

	return {
		step(dt: number) {
			const prev_t = t;
			t += dt;

			const xDot_n = f(prev_t, prev_x);
			const xDot_n_1 = previousXDot[2];
			const xDot_n_2 = previousXDot[1];
			const xDot_n_3 = previousXDot[0];

			const x = sumArrays(
				prev_x,
				scalarMul(
					sumArrays(
						scalarMul(xDot_n, 55),
						scalarMul(xDot_n_1, -59),
						scalarMul(xDot_n_2, 37),
						scalarMul(xDot_n_3, -9)
					),
					(1 / 24) * dt
				)
			);

			previousXDot[0] = xDot_n_2;
			previousXDot[1] = xDot_n_1;
			previousXDot[2] = xDot_n;

			prev_x = x;
			return x;
		},
		reset() {
			t = 0;
			previousXDot = new Array(3).fill(f(t, x0)); // Reset the last 3 derivatives
			initialSolver.reset();
			prev_x = x0;
		}
	};
}
