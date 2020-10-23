//! Module to handle collectives, currently only "value based"
use crate::Pshmem;
use crate::PSHMEM_WORK;
use shmem::collect::Collect;

/// A trait which implement
pub trait ValueCollect<RHS = Self> {
    /// our type
    type T;
    /// reduce and produce maximum value
    fn reduce_max(&self, value: RHS) -> Self::T;
    /// reduce and produce sum
    fn reduce_add(&self, value: RHS) -> Self::T;
}

impl ValueCollect<i64> for Pshmem {
    type T = i64;
    fn reduce_max(&self, value: Self::T) -> Self::T {
        let source = vec![value];
        PSHMEM_WORK.with(|psw| {
            psw.borrow()
                .work_i64
                .shmem_object
                .max_to_all(0, &source)
                .expect("fixme");
            psw.borrow().work_i64.local_part()[0]
        })
    }
    fn reduce_add(&self, value: Self::T) -> Self::T {
        let source = vec![value];
        PSHMEM_WORK.with(|psw| {
            psw.borrow()
                .work_i64
                .shmem_object
                .sum_to_all(0, &source)
                .expect("fixme");
            psw.borrow().work_i64.local_part()[0]
        })
    }
}

impl ValueCollect<usize> for Pshmem {
    type T = usize;
    fn reduce_max(&self, value: Self::T) -> Self::T {
        self.reduce_max(value as i64) as usize
    }
    fn reduce_add(&self, value: Self::T) -> Self::T {
        self.reduce_add(value as i64) as usize
    }
}
