const colorScheme = [
    {
        color: [255, 0, 255],
        members: ['G', 'Y', 'S', 'T', 'N', 'C', 'Q']
    },
    {
        color: [70, 156, 118],
        members: ['V', 'I', 'L', 'P', 'F', 'M', 'W', 'A']
    },
    {
        color: [255, 140, 0],
        members: ['H']
    },
    {
        color: [192, 0, 0],
        members: ['D', 'E']
    },
    {
        color: [0, 0, 255],
        members: ['K', 'R']
    }
]

const unitSize = 10;
const spaceBetweenArchitectures = 250;
// const leftMargin = 60;
const leftMargin = 400;
const rightMargin = 60;
const topMargin = 50;
const canvasWidth = max_length * unitSize + rightMargin + leftMargin;
const canvasHeight = (Object.keys(lengths).length * spaceBetweenArchitectures + topMargin);
const spaceBeforeArchitectureName = 30;
const blackLineThickness = 8;

function setup() {
    background(255);
    createCanvas(canvasWidth, canvasHeight);
    let i = 0;
    for (const architectureName in architectures) {
        let y;
        const familyYOffset = spaceBetweenArchitectures / 2 + topMargin;
        y = familyYOffset + i * spaceBetweenArchitectures;
        // drawArchitectureName(architectureName, y - (76));
        drawArchitectureName(architectureName, y+10);
        drawArchitecture(architectures[architectureName]['architecture_string'], architectures[architectureName]['domain_start_structure'], conservedResidues[architectureName], lengths[architectureName], y);
        i++;
    }
}

function computeXForPositionInDomain(positionInDomain){
    return ((positionInDomain + 0) * unitSize) + leftMargin;
}

function drawArchitecture(architectureString, domain_start_structure, conservedResidues, sequencelength, y) {
    for (const i in architectureString) {
        const positionInDomain = parseInt(i);
        const positionInCompleteProtein = positionInDomain + domain_start_structure;
        const conservedResidue = conservedResidues[positionInCompleteProtein];
        const isConserved = !!conservedResidue;
        const xCoord = computeXForPositionInDomain(positionInDomain)
        switch (architectureString[i]) {
            case 'l':
                drawLoop(xCoord, y, isConserved);
                break;
            case 'h':
                drawHelix(xCoord, y, isConserved);
                break;
            case 's':
                drawSheet(xCoord, y, isConserved);
                break;
            case 'u':
                drawUndefined(xCoord, y, isConserved);
                break;
        }
        if(isConserved){
             drawConservedResidueLetter(conservedResidue, xCoord, y, positionInCompleteProtein)
        }
    }
}

function getConservedResidueColor(conservedResidue) {
    let scheme = colorScheme.filter(({ members }) => members.includes(conservedResidue));
    if (scheme.length === 0) {
        throw "No color scheme defined for " + conservedResidue;
    }
    return scheme[0].color;
}

function getConservedResidueLetter(conservedResidue){
    if(typeof conservedResidue == 'object')
        return conservedResidue.aminoacid;
    return conservedResidue;
}

function getConservedResidueOffset(conservedResidue){
    if(typeof conservedResidue == 'object')
        return conservedResidue.offset;
    return 0;
}

function drawConservedResidueLetter(conservedResidue, xCoord, y, positionInCompleteProtein) {
    const conservedResidueLetter = getConservedResidueLetter(conservedResidue);
    const xOffsetForConservedResidue = getConservedResidueOffset(conservedResidue);

    const xOffsetForPositionNumber = -3 + xOffsetForConservedResidue;
    // const yOffsetForPositionNumber = -50;
    const yOffsetForPositionNumber = 35;
    const xOffsetForLetter = xOffsetForConservedResidue;
    const yOffsetForLetter = -10;
    // const size_residue_text = 50;
    const size_residue_text = 25;
    const color_residue_text = getConservedResidueColor(conservedResidueLetter);
    textSize(size_residue_text);
    fill(...color_residue_text);
    text(conservedResidueLetter, xCoord - (size_residue_text / 4) + xOffsetForLetter, y + yOffsetForLetter);

    const color_position_text = [0, 0, 0];
    const size_position_text = 13;
    textSize(size_position_text);
    fill(...color_position_text);
    text(positionInCompleteProtein, xCoord - size_position_text / 4 + xOffsetForPositionNumber, y - yOffsetForPositionNumber);
}

function drawLoop(x, y, isConserved) {
    fill(200, 200, 200);
    noStroke();
    if (isConserved) {
        fill(0);
    }
    rect(x, y + unitSize * .25, unitSize, unitSize * .5);
}

function drawHelix(x, y, isConserved) {
    fill(130, 130, 130);
    noStroke();
    if (isConserved) {
        fill(0);
    }
    rect(x, y - 6, unitSize, unitSize * 2);
}

function drawSheet(x, y, isConserved) {
    fill(45, 130, 80);
    noStroke();
    if (isConserved) {
        fill(0);
    }
    rect(x, y - 6, unitSize, unitSize * 2)
}

function drawUndefined(x, y, isConserved) {
    fill(200, 200, 200);
    noStroke();
    if (isConserved) {
        fill(0);
    }
    rect(x + unitSize * .25, y + unitSize * .25, unitSize * 0.5, unitSize * .5);
}

function drawArchitectureName(architectureName, y) {
    const size = 30;
    const color = [0, 0, 0];
    textSize(size);
    fill(...color);
    text(architectureName + ':', 0, y);
}

function drawBlackLine(y) {
    const color = [150, 150, 150];
    fill(...color);
    stroke(color)
    rect(0, y, canvasWidth, blackLineThickness)
}

function drawVerticalLine(x) {
    const color = [150, 150, 150];
    fill(...color);
    stroke(color)
    rect(x, blackLineThickness, blackLineThickness, canvasHeight)
}