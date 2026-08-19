import {
	mkdirSync,
	readFileSync,
	renameSync,
	rmSync,
	statSync,
	writeFileSync,
} from "node:fs";
import { homedir } from "node:os";
import { join } from "node:path";
import { randomUUID } from "node:crypto";
import { setTimeout as sleep } from "node:timers/promises";
import type { ExtensionAPI } from "@earendil-works/pi-coding-agent";

type Limit = { requests: number; windowMs: number };
type Reservation = { timestamps: number[]; waitMs: number };

const LOCK_POLL_MS = 25;
const STALE_LOCK_MS = 10_000;

export function planReservation(
	timestamps: number[],
	now: number,
	limit: Limit,
): Reservation {
	const cutoff = now - limit.windowMs;
	const active = timestamps
		.filter((timestamp) => Number.isFinite(timestamp) && timestamp > cutoff && timestamp <= now)
		.sort((left, right) => left - right);

	if (active.length < limit.requests) {
		active.push(now);
		return { timestamps: active, waitMs: 0 };
	}

	const nextExpiry = active[active.length - limit.requests] + limit.windowMs;
	return { timestamps: active, waitMs: Math.max(1, nextExpiry - now) };
}

function loadLimits(): Map<string, Limit> {
	const config = JSON.parse(
		readFileSync(new URL("./aqueduct-rate-limit.json", import.meta.url), "utf8"),
	) as { limits?: Record<string, Partial<Limit>> };
	const limits = new Map<string, Limit>();

	for (const [model, limit] of Object.entries(config.limits ?? {})) {
		if (!Number.isInteger(limit.requests) || (limit.requests ?? 0) < 1)
			throw new Error(`Invalid requests limit for ${model}`);
		if (!Number.isInteger(limit.windowMs) || (limit.windowMs ?? 0) < 1)
			throw new Error(`Invalid windowMs limit for ${model}`);
		limits.set(model, limit as Limit);
	}

	return limits;
}

async function acquireLock(path: string, signal?: AbortSignal): Promise<void> {
	for (;;) {
		signal?.throwIfAborted();
		try {
			mkdirSync(path);
			return;
		} catch (error) {
			if ((error as NodeJS.ErrnoException).code !== "EEXIST") throw error;
		}

		try {
			if (Date.now() - statSync(path).mtimeMs > STALE_LOCK_MS) {
				rmSync(path, { recursive: true, force: true });
				continue;
			}
		} catch (error) {
			if ((error as NodeJS.ErrnoException).code !== "ENOENT") throw error;
		}

		await sleep(LOCK_POLL_MS, undefined, signal ? { signal } : undefined);
	}
}

function readTimestamps(path: string): number[] {
	try {
		const value = JSON.parse(readFileSync(path, "utf8"));
		return Array.isArray(value) ? value.filter((item): item is number => typeof item === "number") : [];
	} catch (error) {
		if ((error as NodeJS.ErrnoException).code === "ENOENT") return [];
		throw error;
	}
}

function writeTimestamps(path: string, timestamps: number[]): void {
	const temporary = `${path}.${process.pid}.${randomUUID()}.tmp`;
	try {
		writeFileSync(temporary, JSON.stringify(timestamps));
		renameSync(temporary, path);
	} finally {
		rmSync(temporary, { force: true });
	}
}

async function reserve(
	stateRoot: string,
	model: string,
	limit: Limit,
	signal?: AbortSignal,
): Promise<void> {
	const bucket = Buffer.from(model).toString("base64url");
	const statePath = join(stateRoot, `${bucket}.json`);
	const lockPath = join(stateRoot, `${bucket}.lock`);

	for (;;) {
		await acquireLock(lockPath, signal);
		let reservation: Reservation;
		try {
			reservation = planReservation(readTimestamps(statePath), Date.now(), limit);
			writeTimestamps(statePath, reservation.timestamps);
		} finally {
			rmSync(lockPath, { recursive: true, force: true });
		}

		if (reservation.waitMs === 0) return;
		await sleep(reservation.waitMs, undefined, signal ? { signal } : undefined);
	}
}

export default function (pi: ExtensionAPI) {
	const limits = loadLimits();
	const agentDir = process.env.PI_CODING_AGENT_DIR?.trim() || join(homedir(), ".pi", "agent");
	const stateRoot = join(agentDir, "cache", "aqueduct-rate-limit");
	mkdirSync(stateRoot, { recursive: true });

	pi.on("before_provider_request", async (_event, ctx) => {
		if (!ctx.model) return;
		const model = `${ctx.model.provider}/${ctx.model.id}`;
		const limit = limits.get(model);
		if (limit) await reserve(stateRoot, model, limit, ctx.signal);
	});
}
