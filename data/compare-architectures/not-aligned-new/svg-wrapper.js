let svgOutput = '';
let fillState = [0,0,0];
let textSizeState = 10;
let textStyleState = "normal";

function fill_svg(...params) {
    if(params.length == 1){
        if(typeof params[0] == "string"){
            fillState[0] = Number("0x" + params[0].substring(1, 3));
            fillState[1] = Number("0x" + params[0].substring(3, 5));
            fillState[2] = Number("0x" + params[0].substring(5, 7));
        }
        else {
            fillState[0] = params[0];
            fillState[1] = params[0];
            fillState[2] = params[0];
        }
    }
    else {
        fillState[0] = params[0];
        fillState[1] = params[1];
        fillState[2] = params[2];
    }
    fill(...params)
}

function rect_svg(...params) {
    let [x, y, width, height] = params;
    svgOutput += `<rect fill="rgb(${fillState[0]}, ${fillState[1]}, ${fillState[2]})" x="${x}" y="${y}" width="${width}" height="${height}" />`
    rect(...params)
}

function background_svg(...params) {
    background(...params)
}

function noStroke_svg(...params) {
    noStroke(...params)
}

function textSize_svg(...params) {
    textSizeState = params[0];
    textSize(...params)
}

function stroke_svg(...params) {
    stroke(...params)
}

function createCanvas_svg(...params) {
    let [width, height] = params;
    svgOutput += `<svg viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg">`;
    createCanvas(...params)
}

function text_svg(...params) {
    let [string, x, y] = params;
    let style = textStyleState == BOLD?'bold':'';
    svgOutput += `<text x="${x}" y="${y}" style="fill: rgb(${fillState[0]}, ${fillState[1]}, ${fillState[2]}); font: ${style} ${textSizeState}px arial;">${string}</text>`;
    text(...params)
}

function textStyle_svg(...params) {
    textStyleState = params[0];
    textStyle(...params)
}