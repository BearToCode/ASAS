export function sumArrays(...arrays: number[][]) {
	return arrays.reduce(
		(accArray, arr) => accArray.map((acc, idx) => acc + arr[idx]),
		Array(arrays[0].length).fill(0)
	);
}

export function scalarMul(a: number[], k: number) {
	return a.map((elem) => elem * k);
}
