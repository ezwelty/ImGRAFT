% SENSORWIDTH  Return the CCD sensor width for a digital camera model.
%
%   width = sensorWidth(makemodel)
%
% Returns the CCD sensor width in mm for the specified camera make and 
% model. Data is from Digital Photography Review (www.dpreview.com) and 
% Noah Snavely's Bundler code (phototour.cs.washington.edu/bundler/).
% 
%
% Input:    makemodel   camera name ([Make ' ' Model])
%           
% Output:   width       sensor width in mm     
%
% See also fmm2fpx, fpx2fmm.

function width = sensorWidth(makemodel)

% Check number of arguments
if (nargin < 1)
    error('Specify make & model of camera.')
end

% Check type of arguments
if (~ischar(makemodel))
    error('Camera make & model must be a character string: [Make Model]');
end

% Camera List
% [deblank(strtrim(Make)) ' ' deblank(strtrim(Model))]
cameras = {
    'Asahi Optical Co.,Ltd. PENTAX Optio330RS';
    'Canon Canon DIGITAL IXUS 400';
    'Canon Canon DIGITAL IXUS 40';
    'Canon Canon DIGITAL IXUS 430';
    'Canon Canon DIGITAL IXUS 500';
    'Canon Canon DIGITAL IXUS 50';
    'Canon Canon DIGITAL IXUS 55';
    'Canon Canon DIGITAL IXUS 60';
    'Canon Canon DIGITAL IXUS 65';
    'Canon Canon DIGITAL IXUS 700';
    'Canon Canon DIGITAL IXUS 750';
    'Canon Canon DIGITAL IXUS 800 IS';
    'Canon Canon DIGITAL IXUS II';
    'Canon Canon EOS 10D';
    'Canon Canon EOS-1D Mark II';
    'Canon Canon EOS-1Ds Mark II';
    'Canon Canon EOS  20D';
    'Canon Canon EOS 20D';
    'Canon Canon EOS 300D DIGITAL';
    'Canon Canon EOS 30D';
    'Canon Canon EOS 350D DIGITAL';
    'Canon Canon EOS 400D DIGITAL';
    'Canon Canon EOS 40D';
    'Canon Canon EOS 5D';
    'Canon Canon EOS 5D Mark II';
    'Canon Canon EOS DIGITAL REBEL';
    'Canon Canon EOS DIGITAL REBEL XT';
    'Canon Canon EOS DIGITAL REBEL XTi';
    'Canon Canon EOS DIGITAL REBEL XS';
    'Canon Canon EOS 1000D';
    'Canon Canon EOS Kiss Digital';
    'Canon Canon IXY DIGITAL 600';
    'Canon Canon PowerShot A10';
    'Canon Canon PowerShot A20';
    'Canon Canon PowerShot A400';
    'Canon Canon PowerShot A40';
    'Canon Canon PowerShot A510';
    'Canon Canon PowerShot A520';
    'Canon Canon PowerShot A530';
    'Canon Canon PowerShot A60';
    'Canon Canon PowerShot A620';
    'Canon Canon PowerShot A630';
    'Canon Canon PowerShot A640';
    'Canon Canon PowerShot A700';
    'Canon Canon PowerShot A70';
    'Canon Canon PowerShot A710 IS';
    'Canon Canon PowerShot A75';
    'Canon Canon PowerShot A80';
    'Canon Canon PowerShot A85';
    'Canon Canon PowerShot A95';
    'Canon Canon PowerShot G1';
    'Canon Canon PowerShot G2';
    'Canon Canon PowerShot G3';
    'Canon Canon PowerShot G5';
    'Canon Canon PowerShot G6';
    'Canon Canon PowerShot G7';
    'Canon Canon PowerShot G9';
    'Canon Canon PowerShot Pro1';
    'Canon Canon PowerShot S110';
    'Canon Canon PowerShot S1 IS';
    'Canon Canon PowerShot S200';
    'Canon Canon PowerShot S2 IS';
    'Canon Canon PowerShot S30';
    'Canon Canon PowerShot S3 IS';
    'Canon Canon PowerShot S400';
    'Canon Canon PowerShot S40';
    'Canon Canon PowerShot S410';
    'Canon Canon PowerShot S45';
    'Canon Canon PowerShot S500';
    'Canon Canon PowerShot S50';
    'Canon Canon PowerShot S60';
    'Canon Canon PowerShot S70';
    'Canon Canon PowerShot S80';
    'Canon Canon PowerShot SD1000';
    'Canon Canon PowerShot SD100';
    'Canon Canon PowerShot SD10';
    'Canon Canon PowerShot SD110';
    'Canon Canon PowerShot SD200';
    'Canon Canon PowerShot SD300';
    'Canon Canon PowerShot SD400';
    'Canon Canon PowerShot SD450';
    'Canon Canon PowerShot SD500';
    'Canon Canon PowerShot SD550';
    'Canon Canon PowerShot SD600';
    'Canon Canon PowerShot SD630';
    'Canon Canon PowerShot SD700 IS';
    'Canon Canon PowerShot SD750';
    'Canon Canon PowerShot SD800 IS';
    'Canon Canon PowerShot SD880 IS';
    'Canon Canon PowerShot SX100 IS';
    'Canon EOS 300D DIGITAL';
    'Canon EOS DIGITAL REBEL';
    'Canon PowerShot A510';
    'Canon PowerShot S30';
    'CASIO COMPUTER CO.,LTD. EX-S500';
    'CASIO COMPUTER CO.,LTD. EX-Z1000';
    'CASIO COMPUTER CO.,LTD  EX-Z30';
    'CASIO COMPUTER CO.,LTD. EX-Z600';
    'CASIO COMPUTER CO.,LTD. EX-Z60';
    'CASIO COMPUTER CO.,LTD  EX-Z750';
    'CASIO COMPUTER CO.,LTD. EX-Z850';
    'CASIO COMPUTER CO.,LTD. EX-Z75';
    'EASTMAN KODAK COMPANY KODAK CX7330 ZOOM DIGITAL CAMERA';
    'EASTMAN KODAK COMPANY KODAK CX7530 ZOOM DIGITAL CAMERA';
    'EASTMAN KODAK COMPANY KODAK DX3900 ZOOM DIGITAL CAMERA';
    'EASTMAN KODAK COMPANY KODAK DX4900 ZOOM DIGITAL CAMERA';
    'EASTMAN KODAK COMPANY KODAK DX6340 ZOOM DIGITAL CAMERA';
    'EASTMAN KODAK COMPANY KODAK DX6490 ZOOM DIGITAL CAMERA';
    'EASTMAN KODAK COMPANY KODAK DX7630 ZOOM DIGITAL CAMERA';
    'EASTMAN KODAK COMPANY KODAK Z650 ZOOM DIGITAL CAMERA';
    'EASTMAN KODAK COMPANY KODAK Z700 ZOOM DIGITAL CAMERA';
    'EASTMAN KODAK COMPANY KODAK Z740 ZOOM DIGITAL CAMERA';
    'FUJIFILM FinePix2600Zoom';
    'FUJIFILM FinePix40i';
    'FUJIFILM FinePix A310';
    'FUJIFILM FinePix A330';
    'FUJIFILM FinePix A600';
    'FUJIFILM FinePix E500';
    'FUJIFILM FinePix E510';
    'FUJIFILM FinePix E550';
    'FUJIFILM FinePix E900';
    'FUJIFILM FinePix F10';
    'FUJIFILM FinePix F30';
    'FUJIFILM FinePix F450';
    'FUJIFILM FinePix F601 ZOOM';
    'FUJIFILM FinePix S3Pro';
    'FUJIFILM FinePix S5000';
    'FUJIFILM FinePix S5200';
    'FUJIFILM FinePix S5500';
    'FUJIFILM FinePix S6500fd';
    'FUJIFILM FinePix S7000';
    'FUJIFILM FinePix Z2';
    'Hewlett-Packard hp 635 Digital Camera';
    'Hewlett-Packard hp PhotoSmart 43x series';
    'Hewlett-Packard HP PhotoSmart 618 (V1.1)';
    'Hewlett-Packard HP PhotoSmart C945 (V01.61)';
    'Hewlett-Packard HP PhotoSmart R707 (V01.00)';
    'KONICA MILOLTA  DYNAX 5D';
    'Konica Minolta Camera, Inc. DiMAGE A2';
    'KONICA MINOLTA CAMERA, Inc. DiMAGE G400';
    'Konica Minolta Camera, Inc. DiMAGE Z2';
    'KONICA MINOLTA DiMAGE A200';
    'KONICA MINOLTA DiMAGE X1';
    'KONICA MINOLTA  DYNAX 5D';
    'Minolta Co., Ltd. DiMAGE F100';
    'Minolta Co., Ltd. DiMAGE Xi';
    'Minolta Co., Ltd. DiMAGE Xt';
    'Minolta Co., Ltd. DiMAGE Z1';
    'NIKON COOLPIX L3';
    'NIKON COOLPIX P2';
    'NIKON COOLPIX S4';
    'NIKON COOLPIX S7c';
    'NIKON COOLPIX S50c';
    'NIKON CORPORATION NIKON D100';
    'NIKON CORPORATION NIKON D1';
    'NIKON CORPORATION NIKON D1H';
    'NIKON CORPORATION NIKON D200';
    'NIKON CORPORATION NIKON D300S';
    'NIKON CORPORATION NIKON D2H';
    'NIKON CORPORATION NIKON D2X';
    'NIKON CORPORATION NIKON D3';
    'NIKON CORPORATION NIKON D3X';
    'NIKON CORPORATION NIKON D40';
    'NIKON CORPORATION NIKON D50';
    'NIKON CORPORATION NIKON D60';
    'NIKON CORPORATION NIKON D70';
    'NIKON CORPORATION NIKON D70s';
    'NIKON CORPORATION NIKON D80';
    'NIKON E2500';
    'NIKON E2500';
    'NIKON E3100';
    'NIKON E3200';
    'NIKON E3700';
    'NIKON E4200';
    'NIKON E4300';
    'NIKON E4500';
    'NIKON E4600';
    'NIKON E5000';
    'NIKON E5200';
    'NIKON E5400';
    'NIKON E5600';
    'NIKON E5700';
    'NIKON E5900';
    'NIKON E7600';
    'NIKON E775';
    'NIKON E7900';
    'NIKON E7900';
    'NIKON E8800';
    'NIKON E990';
    'NIKON E995';
    'NIKON S1';
    'Nokia N80';
    'Nokia N80';
    'Nokia N93';
    'Nokia N95';
    'OLYMPUS CORPORATION C-5000Z';
    'OLYMPUS CORPORATION C5060WZ';
    'OLYMPUS CORPORATION C750UZ';
    'OLYMPUS CORPORATION C765UZ';
    'OLYMPUS CORPORATION C8080WZ';
    'OLYMPUS CORPORATION X250,D560Z,C350Z';
    'OLYMPUS CORPORATION X-3,C-60Z';
    'OLYMPUS CORPORATION X400,D580Z,C460Z';
    'OLYMPUS IMAGING CORP. E-500';
    'OLYMPUS IMAGING CORP. FE115,X715';
    'OLYMPUS IMAGING CORP. SP310';
    'OLYMPUS IMAGING CORP. SP510UZ';
    'OLYMPUS IMAGING CORP. SP550UZ';
    'OLYMPUS IMAGING CORP. uD600,S600';
    'OLYMPUS_IMAGING_CORP. X450,D535Z,C370Z';
    'OLYMPUS IMAGING CORP. X550,D545Z,C480Z';
    'OLYMPUS OPTICAL CO.,LTD C2040Z';
    'OLYMPUS OPTICAL CO.,LTD C211Z';
    'OLYMPUS OPTICAL CO.,LTD C2Z,D520Z,C220Z';
    'OLYMPUS OPTICAL CO.,LTD C3000Z';
    'OLYMPUS OPTICAL CO.,LTD C300Z,D550Z';
    'OLYMPUS OPTICAL CO.,LTD C4100Z,C4000Z';
    'OLYMPUS OPTICAL CO.,LTD C750UZ';
    'OLYMPUS OPTICAL CO.,LTD X-2,C-50Z';
    'OLYMPUS SP550UZ';
    'OLYMPUS X100,D540Z,C310Z';
    'Panasonic DMC-FX01';
    'Panasonic DMC-FX07';
    'Panasonic DMC-FX9';
    'Panasonic DMC-FZ20';
    'Panasonic DMC-FZ2';
    'Panasonic DMC-FZ30';
    'Panasonic DMC-FZ50';
    'Panasonic DMC-FZ5';
    'Panasonic DMC-FZ7';
    'Panasonic DMC-LC1';
    'Panasonic DMC-LC33';
    'Panasonic DMC-LX1';
    'Panasonic DMC-LZ2';
    'Panasonic DMC-TZ1';
    'Panasonic DMC-TZ3';
    'PENTAX Corporation PENTAX *ist DL';
    'PENTAX Corporation PENTAX *ist DS2';
    'PENTAX Corporation PENTAX *ist DS';
    'PENTAX Corporation PENTAX K100D';
    'PENTAX Corporation PENTAX Optio 450';
    'PENTAX Corporation PENTAX Optio 550';
    'PENTAX Corporation PENTAX Optio E10';
    'PENTAX Corporation PENTAX Optio S40';
    'PENTAX Corporation PENTAX Optio S4';
    'PENTAX Corporation PENTAX Optio S50';
    'PENTAX Corporation PENTAX Optio S5i';
    'PENTAX Corporation PENTAX Optio S5z';
    'PENTAX Corporation PENTAX Optio SV';
    'PENTAX Corporation PENTAX Optio WP';
    'RICOH CaplioG3 modelM';
    'RICOH Caplio GX';
    'RICOH Caplio R30';
    'Samsung Digimax 301';
    'Samsung Techwin <Digimax i5, Samsung #1>';
    'SAMSUNG TECHWIN Pro 815';
    'SONY DSC-F828';
    'SONY DSC-N12';
    'SONY DSC-P100';
    'SONY DSC-P10';
    'SONY DSC-P12';
    'SONY DSC-P150';
    'SONY DSC-P200';
    'SONY DSC-P52';
    'SONY DSC-P72';
    'SONY DSC-P73';
    'SONY DSC-P8';
    'SONY DSC-R1';
    'SONY DSC-S40';
    'SONY DSC-S600';
    'SONY DSC-T9';
    'SONY DSC-V1';
    'SONY DSC-W1';
    'SONY DSC-W30';
    'SONY DSC-W50';
    'SONY DSC-W5';
    'SONY DSC-W7';
    'SONY DSC-W80';
    'SONY DCR-HC96';
    'NIKON E8700'
};

% Sensor Widths (mm)
% source: Digital Photography Review
% www.dpreview.com & www.dpreview.com/news/0210/02100402sensorsizes.asp
widths = [
    7.176; 
    7.176;  
    5.76;   
    7.176;  
    7.176;  
    5.76;   
    5.76;   
    5.76;   
    5.76;   
    7.176;  
    7.176;  
    5.76;   
    5.27;   
    22.7;
    28.7;
    35.95;   
    22.5;
    22.5;
    22.66;
    22.5;
    22.2;
    22.2;
    22.2;
    35.8;
    36;
    22.66;
    22.2;
    22.2;
    22.2;
    22.2;
    22.66;
    7.176;  
    5.27;   
    5.27;   
    4.54;   
    5.27;   
    5.76;   
    5.76;   
    5.76;   
    5.27;   
    7.176;  
    7.176;  
    7.176;  
    5.76;   
    5.27;   
    5.76;   
    5.27;   
    7.176;  
    5.27;   
    7.176;  
    7.176;  
    7.176;  
    7.176;  
    7.176;  
    7.176;  
    7.176;  
    7.600;  
    8.8;    
    5.27;   
    5.27;   
    5.27;   
    5.76;   
    7.176;  
    5.76;   
    7.176;  
    7.176;  
    7.176;  
    7.176;  
    7.176;  
    7.176;  
    7.176;  
    7.176;  
    7.176;  
    5.75;   
    5.27;   
    5.75;   
    5.27;   
    5.76;   
    5.76;   
    5.76;   
    5.76;   
    7.176;  
    7.176;  
    5.76;   
    5.76;   
    5.76;   
    5.75;   
    5.76;   
    5.75;
    5.75;
    22.66;
    22.66;
    5.75;   
    7.176;  
    5.76;   
    7.716; 
    5.76;   
    5.76;   
    7.176; 
    7.176; 
    7.176;
    5.75;
    5.27; 
    5.76; 
    7.176; 
    7.176; 
    5.27; 
    5.76; 
    7.176; 
    5.76; 
    5.76; 
    5.76; 
    5.27;   
    7.600;  
    5.27;   
    5.27;   
    7.600;  
    5.76;   
    5.76;   
    7.600;  
    7.78;   
    7.600; 
    7.600;  
    5.76;   
    7.600;  
    23.0;
    5.27;   
    5.76;   
    5.27;   
    7.600;  
    7.600;  
    5.76;   
    4.54; 
    5.27;  
    5.27;  
    7.176; 
    7.176; 
    23.5;
    8.80; 
    5.76; 
    5.76; 
    8.80;   
    7.176;  
    23.5;
    7.176;  
    5.27;   
    5.27;   
    5.27; 
    5.76;   
    7.176;  
    5.76;   
    5.76;   
    5.75;
    23.7;
    23.7;
    23.7;
    23.6;
    23.6;
    23.3;
    23.7;
    36.0;
    35.9;
    23.7;
    23.7;
    23.6;
    23.7;
    23.7;
    23.6;
    5.27;   
    5.27;   
    5.27;   
    5.27;
    5.27;   
    7.176;  
    7.18;
    7.176;  
    5.76;   
    8.80;   
    7.176;  
    7.176;  
    5.76;   
    8.80;   
    7.176;  
    7.176;  
    5.27;   
    7.176;  
    7.176;  
    8.80;   
    7.176;  
    7.176;  
    5.76;   
    5.27;   
    5.27;   
    4.536;  
    5.7;    
    7.176;  
    7.176; 
    5.27;   
    5.76;   
    8.80;   
    5.76; 
    7.176; 
    5.27;  
    17.3;  
    5.76; 
    7.176; 
    5.75;   
    5.76; 
    5.75; 
    5.27; 
    5.76; 
    6.40;  
    5.27;   
    4.54; 
    7.176; 
    5.4;
    7.176;  
    5.27;  
    7.176; 
    5.76;  
    5.27;   
    5.76;   
    5.75;   
    5.76;   
    5.760;  
    4.54;   
    7.176;  
    7.176;  
    5.760;  
    5.76;   
    8.80;   
    5.760;  
    8.50;   
    5.76;   
    5.75;   
    5.68;   
    23.5;
    23.5;
    23.5;
    23.5;
    7.176; 
    7.176; 
    5.76; 
    5.76; 
    5.76; 
    5.76; 
    5.76; 
    5.76; 
    5.76; 
    5.75; 
    5.27;   
    7.176;  
    5.75;   
    5.27;   
    5.76;
    8.80;   
    8.80;   
    7.176;  
    7.176;  
    7.176;  
    7.176;  
    7.176;  
    7.176;   
    5.27;   
    5.27;   
    5.27;
    5.27;   
    21.5;
    5.27;   
    5.760;  
    7.18;
    7.176;  
    7.176;  
    5.760;  
    5.75;   
    7.176;  
    7.176;  
    5.75;   
    4.80;
    8.80;
];

% Check for match
match = find(strcmp(makemodel,cameras));

% If match found, return sensor width
if (~isempty(match))
    width = widths(match);
else
    warning(['No sensor width found for camera "' makemodel '". Returning NaN.']);
    width = NaN;
end
