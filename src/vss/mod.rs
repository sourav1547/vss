pub mod public_parameters;
pub mod keys;
pub mod sigs;
pub mod common;
pub mod transcript;
pub mod ni_vss;
pub mod poly_com_bls;
pub mod poly_com_ed;
pub mod poly_com_mixed;
pub mod poly_com_mixed_ed;
pub mod recon;

pub enum NiVSSError {
    /// Struct to be signed does not serialize correctly.
    DealingError,
}