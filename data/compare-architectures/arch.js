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
        // color: [255, 0, 0],
        color: [255, 140, 0],
        members: ['D', 'E', 'H']
    },
    {
        color: [0, 0, 255],
        members: ['K', 'R']
    }
]

const spaceBetweenArchitectures = 500;
const leftMargin = 150;
const rightMargin = 150;
const topMargin = 100;
const canvasWidth = 5000 + rightMargin + leftMargin;
const canvasHeight = (Object.keys(lengths).length*spaceBetweenArchitectures+topMargin+3-(65));
const spaceBeforeArchitectureName = 60;
const blastLineThickness = 16;
console.log(canvasWidth);

function setup(){
    background(255);
    //noStroke();
    createCanvas(canvasWidth, canvasHeight);
    let i = 0;

    console.log(architectures)
    for(const architectureName in architectures){
        let y;
        const familyYOffset = spaceBetweenArchitectures/2+topMargin;
        y = familyYOffset+i*spaceBetweenArchitectures;
        drawArchitectureName(architectureName, y-(150));
        drawArchitecture(architectures[architectureName], conservedResidues[architectureName], lengths[architectureName], y);
        i++;
    }
    // drawVerticalLine(0)
    // drawVerticalLine(canvasWidth-blastLineThickness)
}

const unitSize = 20;
function drawArchitecture(architectureString, conservedResidues, sequencelength, y){
    for(const i in architectureString){
        if(i==190 && !!conservedResidues[191]){
            console.log('hello');
        }
        const isConserved = !!conservedResidues[parseInt(i)+1];
        const fillConserved = (...v) => !isConserved?fill(...v):fill(0);
        const indent = max_length - sequencelength;
        switch(architectureString[i]){
            case 'l':
                fillConserved(200);
                drawLoop(i*unitSize+leftMargin+indent*unitSize, y);
                break;
            case 'h':
                fillConserved(45, 130, 80);
                drawHelix(i*unitSize+leftMargin+indent*unitSize, y-(10), isConserved);
                break;
            case 's':
                fillConserved(45, 130, 80);
                // drawSheet(i*unitSize+leftMargin+indent*unitSize, y-(10), !!conservedResidues[i]);
                drawSheet(i*unitSize+leftMargin+indent*unitSize, y-(10), isConserved);
                // fillConserved(200);
                // drawLoop(i*unitSize+leftMargin+indent*unitSize, y);
                break;
        }
        drawConservedResidue(conservedResidues[i], i, (i)*unitSize+leftMargin+indent*unitSize, y);
    }
}

function getColor(conservedResidue){
    let scheme = colorScheme.filter(({members})=> members.includes(conservedResidue));
    if(scheme.length === 0){
        throw "No color scheme defined for "+conservedResidue;
    }
    return scheme[0].color;
}

function drawConservedResidue(conservedResidue, position, x, y){
    let conservedResidueAminoacid = conservedResidue;
    let offset = 0;
    if(!conservedResidue){
        return;
    }
    if (typeof conservedResidue == 'object'){
        conservedResidueAminoacid = conservedResidue.aminoacid;
        offset = conservedResidue.offset;
    }
    // const size_residue_text = 180;
    const size_residue_text = 100;
    const color_residue_text = getColor(conservedResidueAminoacid);
    textSize(size_residue_text);
    fill(...color_residue_text);
    text(conservedResidueAminoacid, x-size_residue_text/4+offset-(27), y-(20));

    // const color_position_text = [0, 0, 0];
    // const size_position_text = 100;
    // textSize(size_position_text);
    // fill(...color_position_text);
    // text(position, x-size_residue_text/4-(10)+offset-(27), y-(170));
}

function drawInside(x, y){
    noStroke();
    rect(x, y+unitSize*.5, unitSize, unitSize*.5)
}

function drawOutside(x, y){
    noStroke();
    rect(x, y, unitSize, unitSize*.5);
}

function drawLoop(x,y){
    noStroke();
    rect(x, y+unitSize*.25, unitSize, unitSize*.5);
}

function drawHelix(x, y, hasConservedResidue){
    if(hasConservedResidue){
        rect(x, y, unitSize, unitSize*3);
        return;
    }
    fill(130, 130, 130);
    rect(x, y, unitSize, unitSize*3);
    // let i;
    // for(i = 0; i<unitSize/2; i++){
    //     stroke(map(i, 0, unitSize/2, 0, 255));
    //     line(x+i, y, x+i, y+unitSize*3);
    // }
    // for(; i< unitSize; i++){
    //     stroke(map(i, unitSize/2, unitSize, 255, 0));
    //     line(x+i, y, x+i, y+unitSize*3);
    // }
}

function drawECHelix(x, y, hasConservedResidue){
    if(hasConservedResidue){
        rect(x, y, unitSize, unitSize*3);
        return;
    }
    fill(200, 200, 200);
    rect(x, y, unitSize, unitSize*3);
    // let i;
    // for(i = 0; i<unitSize/2; i++){
    //     stroke(map(i, 0, unitSize/2, 175, 255));
    //     line(x+i, y, x+i, y+unitSize*3);
    // }
    // for(; i< unitSize; i++){
    //     stroke(map(i, unitSize/2, unitSize, 255, 175));
    //     line(x+i, y, x+i, y+unitSize*3);
    // }
}

function drawSheet(x, y){
    const color = [150, 0, 0];
    stroke(color)
    rect(x, y, unitSize, unitSize)
}

function drawFamilyName({familyName, subFamilyName}, y){
    const fontSize = 300;
    fill(0, 0, 0);
    textSize(fontSize);
    const topMargin = 550;
    textStyle(BOLD);
    text(familyName, leftMargin, y-topMargin);
    textSize(fontSize*.4);
    text(subFamilyName, leftMargin+(2.4*fontSize), y-topMargin+(fontSize*.1));
    textStyle(NORMAL);
}

function drawArchitectureName(architectureName, y){
    // const size = 160;
    const size = 100;
    const color = [0, 0, 0];
    textSize(size);
    fill(...color);
    // text(architectureName.substring(0,4), leftMargin+max_length*unitSize+spaceBeforeArchitectureName, y+40);
    // text(architectureName, leftMargin+max_length*unitSize+spaceBeforeArchitectureName, y+40);
    text(architectureName+':', leftMargin, y-40);
}

function drawBlackLine(y){
    const color = [150, 150, 150];
    fill(...color);
    stroke(color)
    rect(0, y, canvasWidth, blastLineThickness)
}

function drawVerticalLine(x){
    const color = [150, 150, 150];
    fill(...color);
    stroke(color)
    rect(x, blastLineThickness, blastLineThickness, canvasHeight)
}