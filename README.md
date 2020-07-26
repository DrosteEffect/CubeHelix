CUBEHELIX Function
==================

CUBEHELIX is the only colormap-generator you will ever need for MATLAB: CUBEHELIX generates colormaps for published or distributed documents as they are very attractive in full color and yet are suitable for grayscale conversion.

CUBEHELIX creates different colormaps using just a few parameters. The standard algorithm offer very attractive colorschemes for online and electronic documents (e.g. PDF), and yet when printed in grayscale they keep exactly the sequence information of the original data. CUBEHELIX function also includes two extra controls over the range and domain of the colormap values, giving a practically unlimited number of colormaps with many different styles: maximally distinct, multi or single hue, and may be suitable for grayscale printing or even simple grayscale.

### Algorithm ###

CUBEHELIX colorschemes consist of nodes along a tapered helix in the RGB color cube, with a continuous increase in perceived intensity (e.g. black->white). Thus the scheme defines attractive colormaps with a huge choice of hue, saturation and brightness, and yet printing a figure (or image) in Black-and-White (e.g. postscript) results in a monotonically increasing grayscale that retains the brightness order of the original colormap. The sequence information of the colormap is retained even in grayscale, which means an attractive colored image can be printed in grayscale and still be informative to the end-user.

The CUBEHELIX algorithm is defined here: http://astron-soc.in/bulletin/11June/289392011.pdf

For more information and examples: http://www.mrao.cam.ac.uk/~dag/CUBEHELIX/

### Examples ###

    % New colors for the COLORMAP example: 
    S = load('spine'); 
    image(S.X)
    colormap(cubehelix) % default parameters
    
    % New colors for the SURF example: 
    [X,Y,Z] = peaks(30); 
    surfc(X,Y,Z) 
    colormap(cubehelix([],0.5,-1.5,1,1,[0.29,0.92])) 
    axis([-3,3,-3,3,-10,5])

### Bonus Function ###

CUBEHELIX_VIEW creates an interactive figure that allows selection of the colorscheme, and that contains two colorbars showing colors of the colormap and the grayscale equivalent.

R2014b or later: CUBEHELIX_VIEW can also update other axes' or figures' colormaps in real time, for example:

    S = load('spine');
    image(S.X)
    cubehelix_view(gca)

### Notes ###

* The original specification (the links above) misnamed the saturation option as "hue". In the CUBEHELIX function the saturation option is named "satn".