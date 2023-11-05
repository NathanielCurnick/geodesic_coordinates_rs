use wasm_bindgen::prelude::*;

pub type Radians = f64;
pub type Degrees = f64;
pub type Metres = f64;
pub type Kilometres = f64;
pub type Minutes = f64;
pub type Seconds = f64;

#[wasm_bindgen]
pub struct LocBearing {
    pub lat: Radians,
    pub lon: Radians,
    pub bearing: Radians,
}

#[wasm_bindgen]
pub struct DistBearing {
    pub distance: Metres,
    pub bearing: Radians,
}
