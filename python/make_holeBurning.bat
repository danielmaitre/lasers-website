python specHoles_alpha.py

(
echo ^<html^>^<head^>^<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/mathjs/3.14.1/math.js"^>^</script^>
echo ^<script type="text/javascript"^>
echo		function gaussian(det, fwhm^) {
echo			return 2/(fwhm^) * math.pow((math.log(2^)/math.pi^),0.5^) * math.exp((-4*math.log(2^)*math.square(det^)^) / math.square(fwhm^)^);
echo		}
echo		function lorentzian(det, lw^) {
echo			return (lw/(2*math.pi^)^) / (math.square(det^) + math.square(lw^)/4^);
echo		}
echo		function maxAndIndex(arr^) {
echo			currMax = arr[0];
echo			maxIdx = 0;
echo			;
echo			for (var i = 1; i ^< arr.length; i++^) {
echo				current = arr[i];
echo				if (current ^> currMax^) {
echo					currMax = current;
echo					maxIdx = i;
echo				}
echo			}
echo			return {maxLz: currMax, maxLzEl: maxIdx};
echo		}
echo ^</script^>
echo ^</head^>^</html^> 
) >> "holeBurning.html"