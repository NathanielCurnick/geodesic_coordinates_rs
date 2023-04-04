pub type Radians = f64;
pub type Degrees = f64;
pub type Metres = f64;
pub type Kilometres = f64;
pub type Minutes = f64;
pub type Seconds = f64;

pub struct LocBearing {
    pub lat: Radians,
    pub lon: Radians,
    pub bearing: Radians,
}

pub struct DistBearing {
    pub distance: Metres,
    pub bearing: Radians,
}
