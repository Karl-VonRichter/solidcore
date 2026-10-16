//! Native oracle for `intersect.curve_class_p5` (`SmCurveClass.cpp`).
//!
//! Flat input: `op shape dim_minus_2 ignore flag1 flag2 flag3 target index scale
//! source_zone angle n_intervals n_solutions optional_xsect_tol`, followed by both
//! classifications' `n+1` boundaries, both classifications' `3*n` point-class codes,
//! then solutions (`type start_t start_other end_t end_other start_gap end_gap`).
//! A negative optional tolerance selects upstream's default. Shapes are line/circle;
//! targets are real line/plane/vertex/face. Counts and native undefined branches are
//! checked before entry. Output contains status, operation outputs, both full interval
//! states and every mutable solution end. No native pointers cross the boundary.

use crate::{call_with_buffer, FfiError};

mod sys {
    extern "C" {
        pub fn sc_intersect_curve_class_p5_run(
            input: *const f64,
            n: i64,
            out: *mut f64,
            cap: i64,
            len: *mut i64,
        ) -> i32;
    }
}

fn invalid() -> FfiError {
    FfiError::InvalidArgument("invalid curve_class_p5 input or undefined native path".into())
}

/// Validates and executes a bounded oracle input.
///
/// Native topology mutation on bare objects, nonempty endpoint fixups, and the
/// null-edge short-range tangency bug are excluded. Valid non-topology insertion,
/// extraction, vertex projection and merge skip paths remain covered.
pub fn run(input: &[f64]) -> Result<Vec<f64>, FfiError> {
    if input.len() < 15 || input.iter().any(|x| !x.is_finite() || x.abs() > 1e21) {
        return Err(invalid());
    }
    let int = |i: usize, max: u32| -> Result<usize, FfiError> {
        let x = input[i];
        if x < 0.0 || x > f64::from(max) || x.fract() != 0.0 {
            return Err(invalid());
        }
        Ok(x as usize)
    };
    let op = int(0, 6)?;
    int(1, 1)?;
    int(2, 1)?;
    for i in 3..7 {
        int(i, 1)?;
    }
    let target = int(7, 3)?;
    int(8, 1)?;
    let n = int(12, 32)?;
    let ns = int(13, 32)?;
    if input.len() != 15 + 2 * (n + 1) + 6 * n + 7 * ns
        || input[9] <= 0.0
        || input[9] > 1e6
        || !(0.0..=1.0).contains(&input[10])
        || !(-1.0..=1.0).contains(&input[14])
    {
        return Err(invalid());
    }
    for side in 0..2 {
        let start = 15 + side * (n + 1);
        let bounds = &input[start..start + n + 1];
        if bounds.iter().any(|x| !(0.0..=1.0).contains(x))
            || bounds.windows(2).any(|w| w[0] >= w[1])
        {
            return Err(invalid());
        }
    }
    let at = 15 + 2 * (n + 1);
    let classes = &input[at..at + 6 * n];
    for &class in classes {
        let valid = class == 100.0
            || class == 110.0
            || (class == 106.0 && target == 0)
            || (class == 107.0 && target == 1)
            || (class == 104.0 && target == 2)
            || (class == 102.0 && target == 3);
        if !valid {
            return Err(invalid());
        }
        if (op == 4 || op == 5) && (101.0..=104.0).contains(&class) {
            return Err(invalid());
        }
    }
    if op == 2 && classes[..3 * n].chunks_exact(3).any(|c| c[1] != 100.0) {
        return Err(invalid());
    }
    if ((op == 0 || op == 1) && n == 0 && ns > 0) || (op == 3 && n == 0) {
        return Err(invalid());
    }
    if op == 6 && input[4] != 0.0 && input[5] != 0.0 {
        return Err(invalid());
    }
    if op == 0 && target > 1 {
        return Err(invalid());
    }
    // Index 1 means target is the source operand and curve is the classified object.
    // Loading a Face at parameter slot 1 would also read its second parameter (slot 2),
    // which this two-variable fixture does not declare: reject it precisely.
    if op == 1 && input[8] == 0.0 && target == 3 {
        return Err(invalid());
    }
    for sol in input[at + 6 * n..].chunks_exact(7) {
        if sol[0] < 0.0
            || sol[0] > 4.0
            || sol[0].fract() != 0.0
            || sol[1..5].iter().any(|x| !(0.0..=1.0).contains(x))
            || sol[5..].iter().any(|x| *x < 0.0 || *x > 1.0)
        {
            return Err(invalid());
        }
        if op == 0
            && input[3] != 0.0
            && sol[0] == 2.0
            && sol[1] <= sol[3]
            && sol[3] - sol[1] < 0.05 * (input[15 + n] - input[15])
        {
            return Err(invalid());
        }
    }
    call_with_buffer(4096, |out, len| {
        // SAFETY: input's exact size, every count/code and all native preconditions
        // above are validated; output and length have caller-owned valid storage.
        unsafe {
            sys::sc_intersect_curve_class_p5_run(
                input.as_ptr(),
                input.len() as i64,
                out.as_mut_ptr(),
                out.len() as i64,
                len,
            )
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejects_malformed_lengths_counts_and_null_edge_tangency() {
        assert!(run(&[]).is_err());
        let mut v = vec![
            0., 0., 1., 1., 0., 0., 0., 0., 0., 1., 0., 0., 1., 1., -1., 0., 1., 0., 1., 100.,
            100., 100., 100., 100., 100., 2., 0., 0., 0.001, 0.001, 0., 0.,
        ];
        assert!(run(&v).is_err());
        v[12] = 1e20;
        assert!(run(&v).is_err());
        v[12] = 1.;
        v.pop();
        assert!(run(&v).is_err());
    }
}
