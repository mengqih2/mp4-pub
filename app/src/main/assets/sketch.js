"use strict";

/*
 * Declare used P5 globals so that standard doesn't complain.
 */

/* global createCanvas, windowWidth, windowHeight */

/* global clear, stroke, background, line, smooth, strokeWeight */

/* global beginShape, endShape, curveVertex, textSize, text */
var LINE_LENGTH = 40;
var TEXT_SIZE = 23;
var ANG_120 = (Math.PI * 2) / 3;
var ANG_90 = Math.PI / 2;
var ANG_60 = ANG_120 / 2;
var ANG_30 = ANG_60 / 2;
var USE_COLOR = true;
/*
 * Entry point for P5.
 */

function setup() {
  // eslint-disable-line no-unused-vars
  createCanvas(windowWidth, windowHeight);
  setupCanvas();
}

function setColorPreference(useColor) {
  //  USE_COLOR = useColor === "true" ? true : false
  USE_COLOR = useColor;
}
/*
 * Accepts a json molecule object and draws it accordingly.
 *
 * Entry point for Android.
 *
 * Expects a json object of the form
    {
      coreCarbons: unsigned,
      substituents: [
        [{ formula: string, bondCount: unsigned, chainLength: unsigned, type: string }],
        ...,
        ]

      ],
      cyclic: bool
    }
 */

function loadMolecule(theData) {
  // eslint-disable-line no-unused-vars
  var molecule = JSON.parse(theData);
  setupCanvas();

  if (molecule.cyclic) {
    drawCyclicMolecule(
      molecule,
      (windowWidth + LINE_LENGTH / 1.5) / 2,
      (windowHeight - LINE_LENGTH / 2) / 2
    );
  } else {
    drawLinearMolecule(
      molecule,
      (windowWidth -
        (molecule.coreCarbons - 1) * LINE_LENGTH * Math.cos(ANG_30)) /
        2,
      (windowHeight + LINE_LENGTH / 1.25) / 2
    );
  }
}
/*
 * Custom function to draw a line on the screen.
 *
 * Returns a point object of the line's end point.
 */

var drawLine = function drawLine(angle, length, x0, y0) {
  smooth();
  strokeWeight(2);
  var deltaX = Math.cos(angle) * length;
  var deltaY = Math.sin(angle) * length;
  var x1 = x0 + deltaX;
  var y1 = y0 + deltaY;
  line(x0, y0, x1, y1);
  smooth();
  beginShape();
  curveVertex(x0, y0);
  curveVertex(x1, y1);
  endShape();
  return {
    x: x1,
    y: y1
  };
};
/*
 * Clear and reset canvas.
 */

var setupCanvas = function setupCanvas() {
  clear();
  stroke(0);
  background(255);
};
/*
 * Custom function to add the chemical formula to the screen.
 */

var drawText = function drawText(textContent, lineAngle, lineLength, x, y) {
  textSize(TEXT_SIZE);
  strokeWeight(0); // Calculations to offset text from line

  var offset = lineLength / 20;
  var deltaX = Math.cos(lineAngle) * offset;
  var deltaY = Math.sin(lineAngle) * offset; // Improve by making this scale with the angle

  if (deltaY > 0) {
    deltaY += (TEXT_SIZE * 3) / 4;
  } else if (deltaY < 0) {
    deltaY -= (TEXT_SIZE * 1) / 4;
  } // Text is left justified instead of center

  if (deltaX < 0) {
    deltaX -= (textContent.length * TEXT_SIZE * 3) / 4; // If perpendicular to x axis
  } else if (Math.abs(deltaX) < 0.001 && text !== "OH") {
    // Can't use deltaX == 0 due to approx vals
    deltaX -= (textContent.length / 3) * TEXT_SIZE;
  }

  text(textContent, x + deltaX, y + deltaY);
};
/*
 * Draws a double bond on the canvas.
 *
 * Returns the end point of the double bond.
 */

var drawDoubleBond = function drawDoubleBond(angle, length, x, y, color) {
  // Sines and cosines necessary to stagger second line properly
  var deltaX = Math.sin(angle) * 6;
  var deltaY = -Math.cos(angle) * 6;
  var newEndPoint = drawLine(angle, length, x, y); // Offset line to be the double bond

  line(x + deltaX, y + deltaY, newEndPoint.x + deltaX, newEndPoint.y + deltaY);

  if (USE_COLOR === true) {
    fill(color);
    stroke(color);
    drawLine(angle + Math.PI, length / 2, newEndPoint.x, newEndPoint.y);
    drawLine(
      angle + Math.PI,
      length / 2,
      newEndPoint.x + deltaX,
      newEndPoint.y + deltaY
    );
  }

  fill(0);
  stroke(0);
  return newEndPoint;
};

var drawAlcohol = function drawAlcohol(angle, lineLength, x, y) {
  var substituentStartingPoint = drawLine(angle, lineLength, x, y);

  if (USE_COLOR === true) {
    fill(255, 0, 0);
    stroke(255, 0, 0);
    drawLine(
      angle + Math.PI,
      lineLength / 2,
      substituentStartingPoint.x,
      substituentStartingPoint.y
    );
  }

  var textContent = Math.cos(angle) < 0 ? "HO" : "OH";
  drawText(
    textContent,
    angle,
    lineLength,
    substituentStartingPoint.x,
    substituentStartingPoint.y
  );
  fill(0);
  stroke(0);
};

var drawHalogen = function drawHalogen(formula, angle, lineLength, x, y) {
  var startingPoint = drawLine(angle, lineLength, x, y);

  if (USE_COLOR === true) {
    if (formula === "F") {
      fill(0, 255, 0);
      stroke(0, 255, 0);
    } else if (formula === "Cl") {
      fill(0, 255, 0);
      stroke(0, 255, 0);
    } else if (formula === "Br") {
      fill(139, 0, 0);
      stroke(139, 0, 0);
    } else if (formula === "He") {
      fill(0, 255, 255);
      stroke(0, 255, 255);
    }

    drawLine(angle + Math.PI, lineLength / 2, startingPoint.x, startingPoint.y);
  }

  drawText(formula, angle, lineLength, startingPoint.x, startingPoint.y);
  fill(0);
  stroke(0);
};

var drawCarbonyl = function drawCarbonyl(angle, lineLength, x, y) {
  var color = "#FF0000";
  var startingPoint = drawDoubleBond(angle, lineLength, x, y, color);

  if (USE_COLOR === true) {
    fill(255, 0, 0);
    stroke(255, 0, 0);
  }

  drawText("O", angle, lineLength, startingPoint.x, startingPoint.y);
  fill(0);
  stroke(0);
};

var drawSubstituents = function drawSubstituents(
  substituents,
  coreBonds,
  angle,
  lineLength,
  x,
  y
) {
  var num = substituents.length;

  for (var i = 0; i < num; i++) {
    var offset = void 0;

    if (num === 1) {
      offset = 0;
    } else if (num === 2) {
      if (coreBonds === 0) {
        offset = 2 * ANG_90;
      } else if (coreBonds === 1) {
        offset = ANG_120;
      } else if (coreBonds === 2) {
        offset = (2 / 3) * ANG_120;

        if (i === 0) {
          angle = angle - offset / 2;
        }
      }
    } else if (num === 3) {
      if (coreBonds === 0) {
        offset = ANG_120;
      } else if (coreBonds === 1) {
        offset = (2 / 3) * ANG_120;

        if (i === 0) {
          angle = angle - offset / 2;
        }
      }
    } else if (num === 4) {
      offset = ANG_90;
    }

    switch (substituents[i].type) {
      case "CARBONYL":
        drawCarbonyl(angle, lineLength, x, y);
        break;

      case "ALCOHOL":
        drawAlcohol(angle, lineLength, x, y);
        break;

      case "HALOGEN":
        drawHalogen(substituents[i].formula, angle, lineLength, x, y);
        break;

      case "ALKANE":
        drawLinearChain(substituents[i].chainLength, angle, lineLength, x, y);
        break;

      default:
        stroke(255, 5, 5); // Denote error

        drawLinearChain(1, angle, lineLength, x, y);
        stroke(0);
        console.log("Unknown substituent type: ".concat(substituents[i].type));
        break;
    }

    angle += offset;
  }
};
/*
 * Draws an alkyl substituent on the canvas with the desired chain length.
 */

var drawLinearChain = function drawLinearChain(
  chainLength,
  angle,
  lineLength,
  x,
  y
) {
  var startingPoint = {
    x: x,
    y: y
  };

  for (var i = 0; i < chainLength; i++) {
    startingPoint = drawLine(
      angle,
      lineLength,
      startingPoint.x,
      startingPoint.y
    );
    angle = i % 2 === 0 ? angle - ANG_60 : angle + ANG_60;
  }

  stroke(0);
};
/*
 * Draws a linear molecule on the screen.
 *
 * Expects a JSON object in the same form described by loadMolecule().
 */

var drawLinearMolecule = function drawLinearMolecule(molecule, x, y) {
  var angle = -ANG_30;
  var startingPoint = {
    x: x,
    y: y
  };

  for (var i = 0; i < molecule.coreCarbons; i++) {
    var substituentAngle = void 0;

    if (i === 0) {
      substituentAngle = ANG_90;
    } else if (i === molecule.coreCarbons - 1) {
      if (angle < 0) {
        substituentAngle = -ANG_30;
      } else {
        substituentAngle = ANG_30 - ANG_120;
      }
    } else if (angle < 0) {
      substituentAngle = ANG_90;
    } else {
      substituentAngle = -ANG_90;
    }

    if (molecule.substituents[i] && molecule.substituents[i].length > 0) {
      drawSubstituents(
        molecule.substituents[i],
        i === 0 || i === molecule.coreCarbons - 1 ? 1 : 2,
        substituentAngle,
        LINE_LENGTH,
        startingPoint.x,
        startingPoint.y
      );
    }

    if (i !== molecule.coreCarbons - 1) {
      startingPoint = drawLine(
        angle,
        LINE_LENGTH,
        startingPoint.x,
        startingPoint.y
      );
      angle = -angle;
    }
  }
};
/*
 * Draws a cyclic molecule on the screen.
 *
 * Expects a JSON object in the same form described by loadMolecule().
 */

var drawCyclicMolecule = function drawCyclicMolecule(molecule, x, y) {
  var angle = ANG_90;
  var startingPoint = {
    x: x,
    y: y
  };
  var delta = (2 * Math.PI) / molecule.coreCarbons;

  for (var i = 0; i < molecule.coreCarbons; i++) {
    var substituentAngle =
      angle + (delta * (molecule.coreCarbons - 2)) / 4 - 2 * ANG_90;

    if (molecule.substituents[i] && molecule.substituents[i].length > 0) {
      drawSubstituents(
        molecule.substituents[i],
        2,
        substituentAngle,
        LINE_LENGTH,
        startingPoint.x,
        startingPoint.y
      );
    }

    startingPoint = drawLine(
      angle,
      LINE_LENGTH,
      startingPoint.x,
      startingPoint.y
    );
    angle += delta;
  }
};
