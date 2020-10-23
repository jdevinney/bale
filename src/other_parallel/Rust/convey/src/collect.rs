//! Conveyor module for collective functions
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Convey, a conveyor library for rust.  For
/// licence information see the LICENSE file in the top level dirctory
/// of the distribution.
use crate::Convey;

#[derive(Debug)]
/// Type of reduction, full or partial
pub enum PType {
    /// 0..num_ranks
    Full,
    /// 0..my_rank
    Partial,
    /// 0..my_rank-1
    Prior,
}

#[derive(Debug)]
/// Initialization value
pub enum IVal {
    /// minimum val for type
    Min,
    /// maximum val for type
    Max,
    /// zero val for type
    Zero,
}

/// Trait for value collectives
pub trait ValueCollect<RHS = Self>
where
    RHS: std::cmp::PartialOrd + std::ops::Add<Output = RHS>,
{
    /// Type of return
    type T;
    /// Generalized reduction, only one that needs to be implemented
    fn reduce(
        &self,
        value: RHS,
        combine_fn: impl Fn(RHS, RHS) -> RHS,
        partial: PType,
        init_value: IVal,
    ) -> Self::T;
    /// Reduce max value
    fn reduce_max(&self, value: RHS) -> Self::T {
        self.reduce(
            value,
            |lhs: RHS, rhs: RHS| if lhs >= rhs { lhs } else { rhs },
            PType::Full,
            IVal::Min,
        )
    }
    /// Reduce min value
    fn reduce_min(&self, value: RHS) -> Self::T {
        self.reduce(
            value,
            |lhs: RHS, rhs: RHS| if lhs <= rhs { lhs } else { rhs },
            PType::Full,
            IVal::Max,
        )
    }
    /// Full summation
    fn reduce_sum(&self, value: RHS) -> Self::T {
        self.reduce(
            value,
            |lhs: RHS, rhs: RHS| lhs.add(rhs),
            PType::Full,
            IVal::Zero,
        )
    }
    /// Partial summation
    fn reduce_partial_sum(&self, value: RHS) -> Self::T {
        self.reduce(
            value,
            |lhs: RHS, rhs: RHS| lhs.add(rhs),
            PType::Partial,
            IVal::Zero,
        )
    }
    /// Prior summation
    fn reduce_prior_sum(&self, value: RHS) -> Self::T {
        self.reduce(
            value,
            |lhs: RHS, rhs: RHS| lhs.add(rhs),
            PType::Prior,
            IVal::Zero,
        )
    }
}

//convenience macro to create definitions by type
macro_rules! decl {
    ($ty: ident) => {
        impl ValueCollect<$ty> for Convey {
            type T = $ty;
            fn reduce(
                &self,
                value: $ty,
                combine_fn: impl Fn($ty, $ty) -> $ty,
                style: PType,
                init: IVal,
            ) -> Self::T {
                let mut ret: $ty = match init {
                    IVal::Min => core::$ty::MIN,
                    IVal::Max => core::$ty::MAX,
                    IVal::Zero => 0 as $ty,
                };
                let start = match style {
                    PType::Full => 0,
                    PType::Partial => self.my_rank,
                    PType::Prior => self.my_rank + 1,
                };
                // special case, if noone has anything to send
                // TODO: maybe make this work
                if self.num_ranks == 1 && start == self.my_rank + 1 {
                    return ret;
                }
                Convey::session()
                    .pull_fn(|item: $ty, _from_rank| {
                        ret = combine_fn(ret, item);
                    })
                    .push_iter((start..self.num_ranks).map(|x| (value, x)))
                    .finish();
                ret
            }
        }
    };
}

decl!(u128);
decl!(u64);
decl!(u32);
decl!(u16);
decl!(u8);
decl!(usize);

decl!(i128);
decl!(i64);
decl!(i32);
decl!(i16);
decl!(i8);
decl!(isize);

decl!(f32);
decl!(f64);

/// Trait for value collectives
pub trait CollectValues<RHS = Self>
where
    RHS: std::cmp::PartialOrd + std::ops::Add<Output = RHS>,
{
    /// Generalized reduction, only one that needs to be implemented
    fn reduce(
        &self,
        combine_fn: impl Fn(RHS, RHS) -> RHS,
        partial: PType,
        finit_value: IVal,
    ) -> RHS;
    /// Reduce max value
    fn reduce_max(&self) -> RHS {
        self.reduce(
            |lhs: RHS, rhs: RHS| if lhs >= rhs { lhs } else { rhs },
            PType::Full,
            IVal::Min,
        )
    }
    /// Reduce min value
    fn reduce_min(&self) -> RHS {
        self.reduce(
            |lhs: RHS, rhs: RHS| if lhs <= rhs { lhs } else { rhs },
            PType::Full,
            IVal::Max,
        )
    }
    /// Full summation
    fn reduce_sum(&self) -> RHS {
        self.reduce(|lhs: RHS, rhs: RHS| lhs.add(rhs), PType::Full, IVal::Zero)
    }
    /// Partial summation
    fn reduce_partial_sum(&self) -> RHS {
        self.reduce(
            |lhs: RHS, rhs: RHS| lhs.add(rhs),
            PType::Partial,
            IVal::Zero,
        )
    }
    /// Prior summation
    fn reduce_prior_sum(&self) -> RHS {
        self.reduce(|lhs: RHS, rhs: RHS| lhs.add(rhs), PType::Prior, IVal::Zero)
    }
}

//convenience macro to create definitions by type
macro_rules! decl1 {
    ($ty: ident) => {
        impl CollectValues<$ty> for $ty {
            fn reduce(
                &self,
                combine_fn: impl Fn($ty, $ty) -> $ty,
                style: PType,
                init: IVal,
            ) -> $ty {
                let mut ret: $ty = match init {
                    IVal::Min => core::$ty::MIN,
                    IVal::Max => core::$ty::MAX,
                    IVal::Zero => 0 as $ty,
                };
                let convey = Convey::new().expect("convey init failed");
                let start = match style {
                    PType::Full => 0,
                    PType::Partial => convey.my_rank,
                    PType::Prior => convey.my_rank + 1,
                };
                // special case, if noone has anything to send
                // TODO: maybe make this work
                if convey.num_ranks == 1 && start == convey.my_rank + 1 {
                    return ret;
                }
                Convey::session()
                    .pull_fn(|item: $ty, _from_rank| {
                        ret = combine_fn(ret, item);
                    })
                    .push_iter((start..convey.num_ranks).map(|x| (*self, x)))
                    .finish();
                ret
            }
        }
    };
}
decl1!(u128);
decl1!(u64);
decl1!(u32);
decl1!(u16);
decl1!(u8);
decl1!(usize);

decl1!(i128);
decl1!(i64);
decl1!(i32);
decl1!(i16);
decl1!(i8);
decl1!(isize);

decl1!(f32);
decl1!(f64);

#[cfg(test)]
mod tests {
    use super::{CollectValues, ValueCollect};
    #[test]
    fn reduce_sum_usize() {
        let mutex = crate::testing_support::TestingMutex::new();
        let sum = mutex.convey.reduce_sum(42_usize);
        assert_eq!(sum, 42 * mutex.convey.num_ranks);
    }
    #[test]
    fn reduce_max_f64() {
        let mutex = crate::testing_support::TestingMutex::new();
        let sum = mutex.convey.reduce_max(42.0);
        assert_eq!(sum, 42.0 * mutex.convey.num_ranks as f64);
    }
    #[test]
    fn reduce_sum_new_usize() {
        let mutex = crate::testing_support::TestingMutex::new();
        let x = 42_usize;
        let sum = x.reduce_sum();
        assert_eq!(sum, 42 * mutex.convey.num_ranks);
    }
}
