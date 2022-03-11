# Plot 2D Spring in Matlab

This project shows how to plot 2D springs in Matlab.

## Tutorial

[`Spring.m`](Spring.m) is a class file which generates coordinates of a 2D spring for plotting. A new spring instance can be created by passing 'spring radius' and 'number of coils'. For example:
```
spr = Spring(5, 10);
```
creates a spring instance `spr` with radius of 5 and coils number of 10. This instance gives coordinates of a spring connecting two arbitrary points by passing the coordinates of the two points into `spr.getSpr(point1, point2)`. For example:
```
[x, y] = spr.getSpr([-1, 5], [4, 2]);
```
returns the coordinates of the spring connecting point (-1, 5) and point (4, 2), and then you could simply visualize the spring by using:
```
plot(x, y);
```
More details can be found in the demo files.

## Demo Files

### [demo_StaticSpring.m](demo_StaticSpring.m)

This demo file plots a simple static spring connecting two arbitrary points.

### [demo_MovingSpring.m](demo_MovingSpring.m)

This demo file makes an animation of three moving points connectted by two different springs.

## Authors

* **Xue-She Wang** - *Initial work* - [Xue-She](https://github.com/Xue-She)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
