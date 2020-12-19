using OpenStreetMapsX

m =  OpenStreetMapX.get_map_data(pth,use_cache = false, trim_to_connected_graph=true);

class2speed = [100.0, 90.0, 90.0, 70.0, 50.0, 40.0, 20.0, 10.0]

open("weights.csv", "w") do io
    println(io, "from,to,w,speed")
    for (i, e) in enumerate(m.e)
        from = m.v[e[1]]
        to = m.v[e[2]]
        w = m.w[from, to]
        println(io, from - 1, ",", to - 1, ",", w, ",", class2speed[m.class[i]])
    end
end

open("nodeloc.csv", "w") do io
    println(io, "id,lat,lon")
    for n in m.n
        lla = LLA(m.nodes[n], m.bounds)
        println(io, m.v[n] - 1, ",", lla.lat, ",", lla.lon)
    end
end
