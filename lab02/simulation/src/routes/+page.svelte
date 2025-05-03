<script lang="ts">
	import { pd } from '$lib/control';
	import { rk4 } from '$lib/rk4';
	import { nonlinearSystem, type SystemParams } from '$lib/system';
	import P5, { type p5, type Sketch } from 'p5-svelte';

	const CartWidth = 230;
	const CartHeight = 50;
	const WheelRadius = 15;
	const WheelOffset = 40;
	const PendulumLength = 250;
	const PendulumWidth = 5;
	const MassWidth = 20;

	const ScaleFactor = 100;
	const DisturbForce = 200;
	const MaxControl = 1000;

	let windowWidth = $state(0);
	let displayPosition = $state(0);

	let controllerActive = $state(true);

	let leftPressed = false;
	let rightPressed = false;

	const Params: SystemParams = {
		M: 3,
		m: 0.1,
		l: 0.75,
		g: 9.8
	};

	const f = nonlinearSystem(Params);

	let x0 = [0, 0, 0, 0];
	let x = $state(x0);

	const Kp_theta = 606;
	const Kd_theta = 47;
	const Kp_x = -79;
	const Kd_x = -46;

	const controller_theta = pd(Kp_theta, Kd_theta);
	const controller_x = pd(Kp_x, Kd_x);

	const controller = (x: number[], thetaRef: number, xRef: number) => {
		const u_theta = controller_theta(thetaRef - x[2], -x[3]);
		const u_x = controller_x(xRef - x[0], -x[1]);
		const u = u_theta + u_x;

		if (u > MaxControl) {
			return MaxControl;
		} else if (u < -MaxControl) {
			return -MaxControl;
		}

		return u;
	};

	const disturb = (t: number) => {
		if (leftPressed) {
			return -DisturbForce;
		} else if (rightPressed) {
			return DisturbForce;
		}
		return 0;
	};

	const thetaRef = 0;
	const xRef = 0;

	const odefun = (t: number, x: number[]) =>
		f(x, controllerActive ? controller(x, thetaRef, xRef) : 0, disturb(t));
	const integrator = rk4(odefun, x0);

	const sketch: Sketch = (p5) => {
		p5.setup = () => {
			p5.createCanvas(window.innerWidth, window.innerHeight);
		};

		p5.draw = () => {
			const dt = p5.deltaTime / 1000;
			// Avoid too large time steps (for example when the window is not focused)
			if (dt > 0.1) {
				return;
			}

			x = integrator.step(dt);

			const theta = x[2];

			function mod(n: number, m: number) {
				return ((n % m) + m) % m;
			}

			displayPosition = mod(x[0] * ScaleFactor + p5.width / 2, p5.width) - p5.width / 2;

			p5.clear();
			p5.translate(p5.width / 2, p5.height / 2);
			// Ground
			p5.line(-p5.width / 2, 0, p5.width / 2, 0);
			drawCart(p5, displayPosition);
			drawPendulum(p5, displayPosition, theta);
		};
	};

	function drawCart(p5: p5, position: number) {
		p5.rectMode(p5.CENTER);
		p5.strokeWeight(2);
		p5.rect(position, -CartHeight / 2 - WheelRadius * 2, CartWidth, CartHeight);

		p5.ellipse(position - CartWidth / 2 + WheelOffset, -WheelRadius, WheelRadius * 2);
		p5.ellipse(position + CartWidth / 2 - WheelOffset, -WheelRadius, WheelRadius * 2);
	}

	function drawPendulum(p5: p5, position: number, theta: number) {
		p5.push();
		p5.translate(position, -CartHeight - WheelRadius * 2);
		p5.rotate(-theta + Math.PI);
		p5.rectMode(p5.CENTER);
		p5.rect(0, PendulumLength / 2, PendulumWidth, PendulumLength);
		p5.ellipse(0, 0, 10);
		p5.rectMode(p5.CORNER);
		p5.ellipse(0, PendulumLength, MassWidth, MassWidth);
		p5.pop();
	}
</script>

<svelte:window
	bind:innerWidth={windowWidth}
	onkeydown={(e) => {
		if (e.key === 'ArrowLeft') {
			leftPressed = true;
		} else if (e.key === 'ArrowRight') {
			rightPressed = true;
		}
	}}
	onkeyup={(e) => {
		if (e.key === 'ArrowLeft') {
			leftPressed = false;
		} else if (e.key === 'ArrowRight') {
			rightPressed = false;
		}
	}}
/>

<button
	class="m-2 cursor-pointer rounded border px-2 py-1"
	onclick={() => {
		integrator.reset();
	}}
>
	Reset
</button>

<button
	onclick={() => (controllerActive = !controllerActive)}
	class="m-2 cursor-pointer rounded border px-2 py-1"
>
	{controllerActive ? 'Deactivate' : 'Activate'} controller
</button>

<main class="h-screen w-screen">
	<P5 parentDivStyle="width: 100%; height: 100%;" {sketch} />
</main>
