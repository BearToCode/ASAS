export function pd(Kp: number, Kd: number) {
	return (e: number, de: number) => e * Kp + de * Kd;
}
