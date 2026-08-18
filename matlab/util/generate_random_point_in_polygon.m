function [lat, lon] = generate_random_point_in_polygon(vertex_lats, vertex_lons)

    % Rejection-sample uniformly from the polygon's bounding box.
    min_lat = min(vertex_lats);
    max_lat = max(vertex_lats);
    min_lon = min(vertex_lons);
    max_lon = max(vertex_lons);

    is_inside = false;

    while ~is_inside

        lat = min_lat + (max_lat - min_lat) * rand;
        lon = min_lon + (max_lon - min_lon) * rand;

        is_inside = inpolygon(lon, lat, vertex_lons, vertex_lats);
    end
end
