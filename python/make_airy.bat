python airy.py

(
echo ^<html^>^<head^>^<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/mathjs/3.14.1/math.js"^>^</script^>
echo ^<script type="text/javascript"^>
echo		function airy(det, finesse^) {
echo			return 1/(1+(4*math.square(finesse^)/(math.square(math.pi^)^)^)*math.square(math.sin(det/2^)^)^);
echo		}
echo ^</script^>
echo ^</head^>^</html^> 
) >> "airy.html"