export function rk4(f: (t: number, x: number[]) => number[], x0: number[]) {
	let t = 0;
	let prev_x = x0;

	return {
		step(dt: number) {
			const prev_t = t;
			t += dt;

			const k_1 = f(prev_t, prev_x);
			const k_2 = f(prev_t + 0.5 * dt, sumArrays(prev_x, scalarMul(k_1, 0.5 * dt)));
			const k_3 = f(prev_t + 0.5 * dt, sumArrays(prev_x, scalarMul(k_2, 0.5 * dt)));
			const k_4 = f(prev_t + dt, sumArrays(prev_x, scalarMul(k_3, dt)));

			const x = sumArrays(prev_x, scalarMul(sumArrays(k_1, k_2, k_3, k_4), (1 / 6) * dt));
			prev_x = x;
			return x;
		},
		reset() {
			t = 0;
			prev_x = x0;
		}
	};
}

function sumArrays(...arrays: number[][]) {
	return arrays.reduce(
		(accArray, arr) => accArray.map((acc, idx) => acc + arr[idx]),
		Array(arrays[0].length).fill(0)
	);
}

function scalarMul(a: number[], k: number) {
	return a.map((elem) => elem * k);
}
