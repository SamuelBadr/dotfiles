import assert from "node:assert/strict";
import { planReservation } from "../extensions/aqueduct-rate-limit.ts";

const limit = { requests: 3, windowMs: 1_000 };

assert.deepEqual(planReservation([9_500, 9_600], 10_000, limit), {
	timestamps: [9_500, 9_600, 10_000],
	waitMs: 0,
});
assert.equal(planReservation([9_500, 9_600, 9_700], 10_000, limit).waitMs, 500);
assert.equal(planReservation([9_100, 9_200, 9_300, 9_400], 10_000, limit).waitMs, 200);
assert.deepEqual(planReservation([9_000, 9_500], 10_000, limit).timestamps, [9_500, 10_000]);

console.log("aqueduct-rate-limit: ok");
