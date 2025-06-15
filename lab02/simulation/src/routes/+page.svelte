<script lang="ts">
	import { ab4 } from '$lib/ab4';
	import { pd } from '$lib/control';
	import { euler } from '$lib/euler';
	import { heun } from '$lib/heun';
	import { rk4 } from '$lib/rk4';
	import { nonlinearSystem, type SystemParams } from '$lib/system';
	import P5, { type p5, type Sketch } from 'p5-svelte';
	import {
		Button,
		Checkbox,
		FpsGraph,
		Pane,
		Point,
		Folder,
		Monitor,
		List,
		type ListOptions,
		ThemeUtils,
		Slider,
		Separator
	} from 'svelte-tweakpane-ui';

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
	let isPaused = $state(false);
	let chosenIntegrator = $state('rk4');
	const integratorOptions: ListOptions<string> = {
		Euler: 'euler',
		Heun: 'heun',
		'Runge-Kutta 4': 'rk4',
		'Adams-Bashforth (AB4)': 'ab4'
	};

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

	let stateX = $derived(x[0]);
	let stateXDot = $derived(x[1]);
	let stateTheta = $derived(x[2]);
	let stateThetaDot = $derived(x[3]);

	const Default_Kp_theta = 606;
	const Default_Kd_theta = 47;
	const Default_Kp_x = -79;
	const Default_Kd_x = -46;

	let gainsX = $state({
		x: Default_Kp_x,
		y: Default_Kd_x
	});
	let gainsTheta = $state({
		x: Default_Kp_theta,
		y: Default_Kd_theta
	});

	let Kp_x = $derived(gainsX.x);
	let Kd_x = $derived(gainsX.y);
	let Kp_theta = $derived(gainsTheta.x);
	let Kd_theta = $derived(gainsTheta.y);

	const controller_theta = $derived(pd(Kp_theta, Kd_theta));
	const controller_x = $derived(pd(Kp_x, Kd_x));

	let u_theta = $state(0);
	let u_x = $state(0);

	const controller = (x: number[], thetaRef: number, xRef: number) => {
		u_theta = controller_theta(thetaRef - x[2], -x[3]);
		u_x = controller_x(xRef - x[0], -x[1]);
		const u = u_theta + u_x;

		if (u > MaxControl) {
			return MaxControl;
		} else if (u < -MaxControl) {
			return -MaxControl;
		}

		return u;
	};

	let stepIntensity = $state(10);
	let disturbanceFrequency = $state(1);
	let disturbanceAmplitude = $state(10);

	let stepDisturbance = $derived(() => stepIntensity);
	let waveDisturbance = $derived((t: number) => {
		return disturbanceAmplitude * Math.sin(disturbanceFrequency * t);
	});

	let selectedDisturbance = $state<'control' | 'wave' | 'step'>('control');
	const disturbanceOptions: ListOptions<'control' | 'wave' | 'step'> = {
		Control: 'control',
		Wave: 'wave',
		Step: 'step'
	};

	let currentDisturbance = $state(0);

	const disturb = (t: number) => {
		if (selectedDisturbance === 'control') {
			if (leftPressed) {
				currentDisturbance = -DisturbForce;
			} else if (rightPressed) {
				currentDisturbance = DisturbForce;
			} else {
				currentDisturbance = 0;
			}
		} else if (selectedDisturbance === 'wave') {
			currentDisturbance = waveDisturbance(t);
		} else if (selectedDisturbance === 'step') {
			currentDisturbance = stepDisturbance();
		}
		return currentDisturbance;
	};

	const thetaRef = 0;
	const xRef = 0;

	const odefun = (t: number, x: number[]) =>
		f(x, controllerActive ? controller(x, thetaRef, xRef) : 0, disturb(t));

	const eulerIntegrator = euler(odefun, x0);
	const heunIntegrator = heun(odefun, x0);
	const rk4Integrator = rk4(odefun, x0);
	const ab4Integrator = ab4(odefun, x0);

	const integrator = $derived.by(() => {
		switch (chosenIntegrator) {
			case 'euler':
				return eulerIntegrator;
			case 'heun':
				return heunIntegrator;
			case 'rk4':
				return rk4Integrator;
			case 'ab4':
				return ab4Integrator;
			default:
				throw new Error(`Unknown integrator: ${chosenIntegrator}`);
		}
	});

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

			if (!isPaused) {
				x = integrator.step(dt);
			}

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

<main class="h-screen w-screen">
	<P5 parentDivStyle="width: 100%; height: 100%;" {sketch} />
</main>

<Pane title="Simulation Controls" theme={ThemeUtils.presets.light}>
	<Folder title="Config">
		<Checkbox bind:value={controllerActive} label="Controller Active" />
		<Button
			on:click={() => {
				integrator.reset();
			}}
			title="Reset"
			label="Reset Simulation"
		/>

		<Point
			bind:value={gainsX}
			expanded={true}
			label="Gains X (Kp, Kd)"
			picker="inline"
			userExpandable={false}
			optionsX={{
				min: -100,
				max: 100,
				step: 1
			}}
			optionsY={{
				min: -100,
				max: 100,
				step: 1
			}}
		/>

		<Button
			on:click={() => {
				gainsX = {
					x: Default_Kp_x,
					y: Default_Kd_x
				};
			}}
			title="Reset"
			label="Reset Gains X"
		/>

		<Point
			bind:value={gainsTheta}
			expanded={true}
			label="Gains Theta (Kp, Kd)"
			picker="inline"
			userExpandable={false}
			optionsX={{
				min: -1000,
				max: 1000,
				step: 1
			}}
			optionsY={{
				min: -1000,
				max: 1000,
				step: 1
			}}
		/>

		<Button
			on:click={() => {
				gainsTheta = {
					x: Default_Kp_theta,
					y: Default_Kd_theta
				};
			}}
			title="Reset"
			label="Reset Gains Theta"
		/>

		<List bind:value={chosenIntegrator} label="Integrator" options={integratorOptions} />
	</Folder>

	<Folder title="Debug">
		<Checkbox bind:value={isPaused} label="Pause Simulation" />
		<FpsGraph interval={50} label="FPS" rows={5} />
	</Folder>

	<Folder title="Disturbances">
		<List bind:value={selectedDisturbance} label="Disturbance Type" options={disturbanceOptions} />
		<Separator />
		<Slider
			bind:value={stepIntensity}
			label="Step Disturbance Intensity"
			min={0}
			max={300}
			step={1}
		/>
		<Separator />
		<Slider
			bind:value={disturbanceFrequency}
			label="Disturbance Frequency"
			min={0.1}
			max={100}
			step={0.1}
		/>
		<Slider
			bind:value={disturbanceAmplitude}
			label="Disturbance Amplitude"
			min={0}
			max={300}
			step={1}
		/>
	</Folder>
</Pane>

<Pane title="System Evolution" theme={ThemeUtils.presets.light}>
	<Folder title="State">
		<Monitor value={stateX} graph label="Cart Position (m)" min={-2} max={+2} />
		<Monitor value={stateX} />
		<Monitor value={stateXDot} graph label="Cart Velocity (m/s)" min={-10} max={10} />
		<Monitor value={stateXDot} />
		<Monitor value={stateTheta} graph label="Pendulum Angle (rad)" min={-3} max={3} />
		<Monitor value={stateTheta} />
		<Monitor
			value={stateThetaDot}
			graph
			label="Pendulum Angular Velocity (rad/s)"
			min={-10}
			max={+10}
		/>
		<Monitor value={stateThetaDot} />
	</Folder>
	<Folder title="Control">
		<Monitor value={u_x} graph label="Controller X" min={-500} max={500} />
		<Monitor value={u_x} />
		<Monitor value={u_theta} graph label="Controller Theta" min={-500} max={500} />
		<Monitor value={u_theta} />
	</Folder>
	<Folder title="Disturbance">
		<Monitor
			value={currentDisturbance}
			interval={selectedDisturbance === 'control' ? 0.1 : 0}
			graph
			label="Disturbance Force"
			min={-300}
			max={300}
		/>
		<Monitor value={currentDisturbance} />
	</Folder>
</Pane>
