///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Convey, a conveyor library for rust.  For
/// licence information see the LICENSE file in the top level dirctory
/// of the distribution.
extern crate bindgen;

use std::env;
use std::path::PathBuf;

fn main() {
    // Tell cargo to tell rustc to link the system bzip2
    // shared library.
    println!("cargo:rustc-link-lib=sma");

    // Tell cargo to invalidate the built crate whenever either the
    // wrapper changes or we change which shmem we are using
    println!("cargo:rerun-if-changed=src/wrapper.h");
    println!("cargo:rerun-if-env-changed=SHMEM_PATH");

    let shmem_path = match env::var("CRAY_SHMEM_DIR") {
        Ok(e) => {
            println!("cargo:rustc-cfg=cray");
            e
        }
        Err(_) => {
            let path = match env::var("SHMEM_PATH") {
                Ok(e) => e,
                Err(_) => "/usr/local".to_string(),
            };
            println!("cargo:rustc-link-lib=pmi_simple");
            println!("cargo:rustc-link-search={}/lib", path);
            path
        }
    };

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // Get rid of longdoubles, etc. they are not supported by clean ABI
        .blacklist_function("*shmem_longdouble.*")
        .blacklist_function("*shmem_ctx_longdouble.*")
        .blacklist_function("*shmem_float128.*")
        .blacklist_function("*shmem_ctx_float128.*")
        .blacklist_function("*shmem_ld80.*")
        .blacklist_function("*shmem_ctx_ld80.*")
        // The input header we would like to generate
        // bindings for.
        .header("src/wrapper.h")
        .clang_arg("-I")
        .clang_arg(format!("{}/include", shmem_path))
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        //.parse_callbacks(Box::new(bindgen::CargoCallbacks))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
