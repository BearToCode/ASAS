export function heun(f: (t: number, x: number[]) => number[], x0: number[]) {
	let t = 0;
	let prev_x = x0;

	return {
		step(dt: number) {
			const prev_t = t;
			t += dt;

			const k_1 = f(prev_t, prev_x);
			const x_temp = prev_x.map((xi, i) => xi + k_1[i] * dt);
			const k_2 = f(prev_t + dt, x_temp);

			const x = prev_x.map((xi, i) => xi + (k_1[i] + k_2[i]) * (dt / 2));
			prev_x = x;
			return x;
		},
		reset() {
			t = 0;
			prev_x = x0;
		}
	};
}
