import { sumArrays, scalarMul } from './utils';

export function euler(f: (t: number, x: number[]) => number[], x0: number[]) {
	let t = 0;
	let prev_x = x0;

	return {
		step(dt: number) {
			const prev_t = t;
			t += dt;

			const k = f(prev_t, prev_x);
			const x = sumArrays(prev_x, scalarMul(k, dt));
			prev_x = x;
			return x;
		},
		reset() {
			t = 0;
			prev_x = x0;
		}
	};
}
