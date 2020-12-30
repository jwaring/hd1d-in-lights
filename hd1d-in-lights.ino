/*
 * Copyright (C) 2020, SoftWaring Solutions ATF The Miss Trust
 */

#include <math.h>
#include <FastLED.h>


#define G (9.81)               // Gravity.
#define NUM_LEDS 150
#define NCELLS 75
#define LLENGTH (32)
#define UNDEFINED -999
#define DELAY_PIN A0
#define TOGGLE_PIN 2
#define WAVE_PIN 3
#define DELAY_PIN 4
#define BRIGHTNESS_PIN 5
#define LED_PIN 7

// Static variables
uint8_t lambdaLength = LLENGTH;
uint8_t nx = NCELLS;           // Length of channel.
float offset = 0.0;           // Start location of grid.
float dx = 125.0;             // Grid width.
float dt = 1.0;               // Time step in s.
float h = 150.0;              // Depth.
float fCoeff = 0.000;        // Friction coefficient.
float lscale = ((float)NCELLS) / ((float)NUM_LEDS);

// Dynamic variables
float u[NCELLS + 1];    // Current.
float eta[NCELLS + 1];  // Surface elevation.
boolean flags[NCELLS + 1]; // Masking flags.
float cscale, cupper, clower;
uint8_t wave_mode = 0;
uint8_t duration_mode = 0;
uint8_t brightness_mode = 0;
int duration;
float brightness;
float t;                         // Current time.
long stepcount;
int nxp1, nxm1;                 // Number of times step was called.
float (*fn_etaspec)(long step, int i);

CRGB leds[NUM_LEDS];
CRGBPalette16 currentPalette;
TBlendType currentBlending;

/*
   A 1-d Hydrodyanmic model.

   This is a channel model, which implements depth-averaged momentum and continuity equations.
   The channel is represented as an equispaced grid (dx). The surface evelation is forced at either
   end of the channel, using the etaspec function (see config). Hydrodynamic models are time-stepping,
   but the delta time (dt in this case) must not be less than the time it take a parcel of water moving
   at the maximum velocity to traverse a cell dimension (dx). To select an appropriate dt, consider
   the courant number (see https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition)

   @author Jason Waring
*/
float hd_cosine_eta_spec(long step, int i) {
  if (step < 2 * lambdaLength) {
    return 0.5 * (1.0 + cos((M_2_PI * step) / lambdaLength + M_PI));
  } else {
    return UNDEFINED;
  }
}


float hd_triangle_eta_spec(long step, int i) {
  if (step < 2 * lambdaLength) {
    /* Triangular forcing function */
    float ll2 = lambdaLength / 2.0;
    if (step >= ll2) {
      return (lambdaLength - step) / ll2;
    } else {
      return step / ll2;
    }
  } else {
    return UNDEFINED;
  }
}

float hd_square_eta_spec(long step, int i) {
  if (step < 2 * lambdaLength) {
    return 1.0;
  } else {
    return UNDEFINED;
  }
}

float hd_phased_wave_eta_spec (long step, int i) {
  if (step < 2 * lambdaLength) {
    if (i == 0) {
      /* Cosine forcing function of three 'overlaid' waves */
      return 0.5 * (1.0 + cos((M_PI * step) / lambdaLength + M_PI));
    } else {
      return 0.75 * (1.0 + cos((M_2_PI * step) / lambdaLength + M_PI));
    }
  } else {
    return UNDEFINED;
  }
};

void set_wave_mode(uint8_t mode) {
  wave_mode = mode % 4;

  switch (wave_mode) {
    case 0:
      fn_etaspec = hd_square_eta_spec;
      cupper = 2.0;
      clower = -1.0;
      break;

    case 1:
      fn_etaspec = hd_cosine_eta_spec;
      cupper = 0.5;
      clower = -0.5;
      break;

    case 2:
      fn_etaspec = hd_triangle_eta_spec;
      cupper = 1.5;
      clower = -3.0;
      break;

    case 3:
    default:
      fn_etaspec = hd_phased_wave_eta_spec;
      cupper = 1.5;
      clower = -1.0;
      break;
  }

  cscale = 255.0 / (cupper - clower);
}


void set_duration(uint8_t mode) {
  duration_mode = mode % 5;
  switch (duration_mode) {
    case 0:
      duration = 1;
      break;

    case 1:
      duration = 10;
      break;

    case 2:
      duration = 50;
      break;

    case 3:
      duration = 100;
      break;

    case 4:
      duration = 200;
      break;
  }
}

void set_brightness(uint8_t mode) {
  brightness_mode = mode % 5;
  switch (brightness_mode) {
    case 0:
      brightness = 1;
      break;

    case 1:
      brightness = 0.75;
      break;

    case 2:
      brightness = 0.5;
      break;

    case 3:
      brightness = 0.25;
      break;

    case 4:
      brightness = 0.05;
      break;
  }
}


void hd_init(int mode) {
  t = 0;                    // Current time.
  stepcount = 0;            // Number of times step was called.
  nxp1 = nx + 1;            // Length plus one.
  nxm1 = nx - 1;            // Length minus one.

  set_wave_mode(mode);

  // Initialise all arrays
  for (int i = 0; i < nxp1; ++i) {
    u[i] = 0.0;
    eta[i] = 0.0;
    flags[i] = false;
  }

  // Should be set externally, but this is a channel after all.
  flags[0] = true;
  flags[nxm1] = true;
}

// Step forward by one time-step.
void step() {
  // Compute u
  for (int i = 1; i < nx; ++i) {
    float dudt = -(G * (eta[i] - eta[i - 1]) / dx);
    dudt = dudt - u[i] * fCoeff;            // Linear friction.
    u[i] += dt * dudt;
  }

  // U Boundary
  // u[0] = 0.0;
  // u[nx] = 0.0;
  u[0] = u[1];
  u[nx] = 0.0;
  //        u[nx] = u[nxm1];


  // Compute eta
  for (int i = 0; i < nx; ++i) {
    float detadt = -(h * (u[i + 1] - u[i]) / dx);
    eta[i] += dt * detadt;
  }

  // ETA Boundary

  uspec();
  etaspec();

  ++stepcount;
  t += dt;
}

// Manages the U specification points.
void uspec() {
}

// Manages the ETA specification points.
void etaspec() {
  for (int i = 0; i < nx; ++i) {
    if (flags[i]) {
      float e = fn_etaspec(stepcount, i);
      if (e != UNDEFINED) {
        eta[i] = e;
      }
    }
  }
}

void clear_leds() {
  for (int i = 0; i < NUM_LEDS; ++i) {
    leds[i] = CRGB::Black;
  }

  FastLED.show();
}



void setup() {
  hd_init(0);
  set_duration(0);;
  set_brightness(0);

  pinMode(TOGGLE_PIN, INPUT_PULLUP);
  pinMode(WAVE_PIN, INPUT_PULLUP);
  pinMode(DELAY_PIN, INPUT_PULLUP);
  pinMode(BRIGHTNESS_PIN, INPUT_PULLUP);

  FastLED.addLeds<WS2812, LED_PIN, GRB>(leds, NUM_LEDS);
  currentPalette = RainbowColors_p;
  currentBlending = LINEARBLEND;
//  Serial.begin(115200);

  clear_leds();

  delay(2000);
}

void loop() {

  int pinValue = digitalRead(TOGGLE_PIN);
//  Serial.println(pinValue);
  if (!pinValue) {

    if (!digitalRead(WAVE_PIN)) {
      clear_leds();
      delay(1000);  // debounce
      hd_init(wave_mode + 1);

    } else if (!digitalRead(DELAY_PIN)) {
      set_duration(duration_mode + 1);

    } else if (!digitalRead(BRIGHTNESS_PIN)) {
      set_brightness(brightness_mode + 1);
    }
  } else {
    step();
  }

  for (int i = 0; i < NUM_LEDS; ++i) {
    float fi = i * lscale;
    int il = (int)fi;
    int iu = (il >= NUM_LEDS) ? il : il + 1;
    float f = (fi - il);
    float value = eta[il] * (1 - f) + eta[iu] * f;
    float index = max(min(cscale * (value - clower), 255), 0);
    uint8_t ci =  (uint8_t)255 - index;
    leds[i] = ColorFromPalette(currentPalette, ci, (uint8_t)(brightness * max(ci, 1)), currentBlending);
  }

  FastLED.show();

  delay(duration);
}
