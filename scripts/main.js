var map;
var marker;
var poly;

$(document).ready(function () {
	$('#path-checkbox').change(function () {
		if (this.checked) {
			// Turn on ground track
			poly.setMap(map);
		} else {
			// Turn off
			poly.setMap(null);
		}
	});

});

function initialize() {
	var mapOptions = {
		center : {
			lat : 0,
			lng : 0
		},
		zoom : 2
	};

	map = new google.maps.Map(document.getElementById('map-canvas'),
			mapOptions);

	var LatLng = new google.maps.LatLng(0, 0);
	marker = new google.maps.Marker({
			position : LatLng,
			map : map,
			title : 'Satellite'
		});

	polyOptions = {
		strokeColor : '#000000',
		strokeOpacity : 1.0,
		strokeWeight : 3
	}

	poly = new google.maps.Polyline(polyOptions);
	poly.setMap(map);
}

function updateMarker() {
	$.post('main.py', {}, function (json) {
		var parsed = JSON.parse(json);

		var LatLng = new google.maps.LatLng(parsed.latitude, parsed.longitude);
		marker.setPosition(LatLng);

		// Extend the polyline
		var path = poly.getPath();
		path.push(LatLng);

		if ($('#centre-frame').is(":checked")) {
			// Set map centre again
			map.setCenter(LatLng);
		}
	});

}

google.maps.event.addDomListener(window, 'load', initialize);
setInterval(updateMarker, 10000); // Every 10 seconds we should update position