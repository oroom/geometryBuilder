/*:
 # Hello SceneKit
 
 A quick demo of SceneKit, including the new physics-based rendering.
 
 Requires either:
 
 + iPad Playgrounds iOS 10 beta 4+ 
 + or Xcode beta 8
 
 NB when running in Xcode, you need to explicitly set the size of the view (see lines 104-105). Also PBR (physics based rendering) material properties don't show up (because Xcode iOS Simulator isn't backed by Metal?)
 
 */

import SpriteKit // we're just using SpriteKit to create a quick noise texture
import PlaygroundSupport // needed to create the live view
import Foundation
import SceneKit

extension SCNVector3
{
    /**
     * Negates the vector described by SCNVector3 and returns
     * the result as a new SCNVector3.
     */
    func negate() -> SCNVector3 {
        return self * -1
    }
    
    /**
     * Negates the vector described by SCNVector3
     */
    mutating func negated() -> SCNVector3 {
        self = negate()
        return self
    }
    
    /**
     * Returns the length (magnitude) of the vector described by the SCNVector3
     */
    func length() -> Float {
        return sqrtf(x*x + y*y + z*z)
    }
    
    /**
     * Returns the length (magnitude) of the vector described by the SCNVector3
     */
    func lengthSquared() -> Float {
        return x*x + y*y + z*z
    }
    
    /**
     * Normalizes the vector described by the SCNVector3 to length 1.0 and returns
     * the result as a new SCNVector3.
     */
    func normalized() -> SCNVector3 {
        return self / length()
    }
    
    /**
     * Normalizes the vector described by the SCNVector3 to length 1.0.
     */
    mutating func normalize() -> SCNVector3 {
        self = normalized()
        return self
    }
    
    /**
     * Calculates the distance between two SCNVector3. Pythagoras!
     */
    func distance(_ vector: SCNVector3) -> Float {
        return (self - vector).length()
    }
    
    /**
     * Calculates the distance between two SCNVector3. Pythagoras!
     */
    func distanceSquared(_ vector: SCNVector3) -> Float {
        return (self - vector).lengthSquared()
    }
    
    /**
     * Calculates the dot product between two SCNVector3.
     */
    func dot(_ vector: SCNVector3) -> Float {
        return x * vector.x + y * vector.y + z * vector.z
    }
    
    /**
     * Calculates the cross product between two SCNVector3.
     */
    func cross(_ vector: SCNVector3) -> SCNVector3 {
        return SCNVector3Make(y * vector.z - z * vector.y, z * vector.x - x * vector.z, x * vector.y - y * vector.x)
    }
}

/**
 * Adds two SCNVector3 vectors and returns the result as a new SCNVector3.
 */
func + (left: SCNVector3, right: SCNVector3) -> SCNVector3 {
    return SCNVector3Make(left.x + right.x, left.y + right.y, left.z + right.z)
}

/**
 * Increments a SCNVector3 with the value of another.
 */
func += (left: inout SCNVector3, right: SCNVector3) {
    left = left + right
}

func == (left:SCNVector3, right:SCNVector3) -> Bool {
    return left.x == right.x && left.y == right.y && left.z == right.z
}

func != (left:SCNVector3, right:SCNVector3) -> Bool {
    return left.x != right.x || left.y != right.y || left.z != right.z
}

/**
 * Subtracts two SCNVector3 vectors and returns the result as a new SCNVector3.
 */
func - (left: SCNVector3, right: SCNVector3) -> SCNVector3 {
    return SCNVector3Make(left.x - right.x, left.y - right.y, left.z - right.z)
}

/**
 * Decrements a SCNVector3 with the value of another.
 */
func -= (left: inout SCNVector3, right: SCNVector3) {
    left = left - right
}

/**
 * Multiplies two SCNVector3 vectors and returns the result as a new SCNVector3.
 */
func * (left: SCNVector3, right: SCNVector3) -> SCNVector3 {
    return SCNVector3Make(left.x * right.x, left.y * right.y, left.z * right.z)
}

/**
 * Multiplies a SCNVector3 with another.
 */
func *= (left: inout SCNVector3, right: SCNVector3) {
    left = left * right
}

/**
 * Multiplies the x, y and z fields of a SCNVector3 with the same scalar value and
 * returns the result as a new SCNVector3.
 */
func * (vector: SCNVector3, scalar: Float) -> SCNVector3 {
    return SCNVector3Make(vector.x * scalar, vector.y * scalar, vector.z * scalar)
}

/**
 * Multiplies the x and y fields of a SCNVector3 with the same scalar value.
 */
func *= (vector: inout SCNVector3, scalar: Float) {
    vector = vector * scalar
}

/**
 * Divides two SCNVector3 vectors abd returns the result as a new SCNVector3
 */
func / (left: SCNVector3, right: SCNVector3) -> SCNVector3 {
    return SCNVector3Make(left.x / right.x, left.y / right.y, left.z / right.z)
}

/**
 * Divides a SCNVector3 by another.
 */
func /= (left: inout SCNVector3, right: SCNVector3) {
    left = left / right
}

/**
 * Divides the x, y and z fields of a SCNVector3 by the same scalar value and
 * returns the result as a new SCNVector3.
 */
func / (vector: SCNVector3, scalar: Float) -> SCNVector3 {
    return SCNVector3Make(vector.x / scalar, vector.y / scalar, vector.z / scalar)
}

/**
 * Divides the x, y and z of a SCNVector3 by the same scalar value.
 */
func /= (vector: inout SCNVector3, scalar: Float) {
    vector = vector / scalar
}

/**
 * Negate a vector
 */
func SCNVector3Negate(_ vector: SCNVector3) -> SCNVector3 {
    return vector * -1
}

/**
 * Returns the length (magnitude) of the vector described by the SCNVector3
 */
func SCNVector3Length(_ vector: SCNVector3) -> Float
{
    return sqrtf(vector.x*vector.x + vector.y*vector.y + vector.z*vector.z)
}

/**
 * Returns the length squared (magnitude) of the vector described by the SCNVector3
 */
func SCNVector3LengthSquared(_ vector: SCNVector3) -> Float
{
    return (vector.x*vector.x + vector.y*vector.y + vector.z*vector.z)
}

/**
 * Returns the distance between two SCNVector3 vectors
 */
func SCNVector3Distance(_ vectorStart: SCNVector3, vectorEnd: SCNVector3) -> Float {
    return SCNVector3Length(vectorEnd - vectorStart)
}

/**
 * Returns the distance between two SCNVector3 vectors
 */
func SCNVector3Normalize(_ vector: SCNVector3) -> SCNVector3 {
    return vector / SCNVector3Length(vector)
}

/**
 * Calculates the dot product between two SCNVector3 vectors
 */
func SCNVector3DotProduct(_ left: SCNVector3, right: SCNVector3) -> Float {
    return left.x * right.x + left.y * right.y + left.z * right.z
}

/**
 * Calculates the cross product between two SCNVector3 vectors
 */
func SCNVector3CrossProduct(_ left: SCNVector3, right: SCNVector3) -> SCNVector3 {
    return SCNVector3Make(left.y * right.z - left.z * right.y, left.z * right.x - left.x * right.z, left.x * right.y - left.y * right.x)
}

/**
 * Calculates the SCNVector from lerping between two SCNVector3 vectors
 */
func SCNVector3Lerp(_ vectorStart: SCNVector3, vectorEnd: SCNVector3, t: Float) -> SCNVector3 {
    return SCNVector3Make(vectorStart.x + ((vectorEnd.x - vectorStart.x) * t), vectorStart.y + ((vectorEnd.y - vectorStart.y) * t), vectorStart.z + ((vectorEnd.z - vectorStart.z) * t))
}

/**
 * Project the vector, vectorToProject, onto the vector, projectionVector.
 */
func SCNVector3Project(_ vectorToProject: SCNVector3, projectionVector: SCNVector3) -> SCNVector3 {
    let scale: Float = SCNVector3DotProduct(projectionVector, right: vectorToProject) / SCNVector3DotProduct(projectionVector, right: projectionVector)
    let v: SCNVector3 = projectionVector * scale
    return v
}

class GeometryBuildTools {
    
    static func buildBox() -> SCNGeometry {
        let halfSideX:Float = 0.5
        let halfSideY:Float = 0.5
        let halfSideZ:Float = 2.0
        
        let positions = [
            SCNVector3Make(-halfSideX, -halfSideY,  halfSideZ),
            SCNVector3Make( halfSideX, -halfSideY,  halfSideZ),
            SCNVector3Make(-halfSideX, -halfSideY, -halfSideZ),
            SCNVector3Make( halfSideX, -halfSideY, -halfSideZ),
            SCNVector3Make(-halfSideX,  halfSideY,  halfSideZ),
            SCNVector3Make( halfSideX,  halfSideY,  halfSideZ),
            SCNVector3Make(-halfSideX,  halfSideY, -halfSideZ),
            SCNVector3Make( halfSideX,  halfSideY, -halfSideZ),
            
            // repeat exactly the same
            SCNVector3Make(-halfSideX, -halfSideY,  halfSideZ),
            SCNVector3Make( halfSideX, -halfSideY,  halfSideZ),
            SCNVector3Make(-halfSideX, -halfSideY, -halfSideZ),
            SCNVector3Make( halfSideX, -halfSideY, -halfSideZ),
            SCNVector3Make(-halfSideX,  halfSideY,  halfSideZ),
            SCNVector3Make( halfSideX,  halfSideY,  halfSideZ),
            SCNVector3Make(-halfSideX,  halfSideY, -halfSideZ),
            SCNVector3Make( halfSideX,  halfSideY, -halfSideZ),
            
            // repeat exactly the same
            SCNVector3Make(-halfSideX, -halfSideY,  halfSideZ),
            SCNVector3Make( halfSideX, -halfSideY,  halfSideZ),
            SCNVector3Make(-halfSideX, -halfSideY, -halfSideZ),
            SCNVector3Make( halfSideX, -halfSideY, -halfSideZ),
            SCNVector3Make(-halfSideX,  halfSideY,  halfSideZ),
            SCNVector3Make( halfSideX,  halfSideY,  halfSideZ),
            SCNVector3Make(-halfSideX,  halfSideY, -halfSideZ),
            SCNVector3Make( halfSideX,  halfSideY, -halfSideZ)
        ]
        
        let normals = [
            SCNVector3Make( 0, -1, 0),
            SCNVector3Make( 0, -1, 0),
            SCNVector3Make( 0, -1, 0),
            SCNVector3Make( 0, -1, 0),
            
            SCNVector3Make( 0, 1, 0),
            SCNVector3Make( 0, 1, 0),
            SCNVector3Make( 0, 1, 0),
            SCNVector3Make( 0, 1, 0),
            
            
            SCNVector3Make( 0, 0,  1),
            SCNVector3Make( 0, 0,  1),
            SCNVector3Make( 0, 0, -1),
            SCNVector3Make( 0, 0, -1),
            
            SCNVector3Make( 0, 0, 1),
            SCNVector3Make( 0, 0, 1),
            SCNVector3Make( 0, 0, -1),
            SCNVector3Make( 0, 0, -1),
            
            
            SCNVector3Make(-1, 0, 0),
            SCNVector3Make( 1, 0, 0),
            SCNVector3Make(-1, 0, 0),
            SCNVector3Make( 1, 0, 0),
            
            SCNVector3Make(-1, 0, 0),
            SCNVector3Make( 1, 0, 0),
            SCNVector3Make(-1, 0, 0),
            SCNVector3Make( 1, 0, 0),
            ]
        
        let indices:[CInt] = [
            // bottom
            0, 2, 1,
            1, 2, 3,
            // back
            10, 14, 11,  // 2, 6, 3,   + 8
            11, 14, 15,  // 3, 6, 7,   + 8
            // left
            16, 20, 18,  // 0, 4, 2,   + 16
            18, 20, 22,  // 2, 4, 6,   + 16
            // right
            17, 19, 21,  // 1, 3, 5,   + 16
            19, 23, 21,  // 3, 7, 5,   + 16
            // front
            8,  9, 12,  // 0, 1, 4,   + 8
            9, 13, 12,  // 1, 5, 4,   + 8
            // top
            4, 5, 6,
            5, 7, 6
        ]
        
        let colors = [vector_float4](repeating: vector_float4(1.0,1.0,1.0,0.0), count: positions.count)
        
        let geom = packGeometryTriangles(positions, indexList: indices, colourList: colors, normalsList: normals)
        return geom
    }
    
    static func buildGeometryWithColor(_ geometry:SCNGeometry, color:vector_float4) -> SCNGeometry {
        
        let indexElement = geometry.element(at: 0)
        
        let vertexSource = geometry.sources(for: SCNGeometrySource.Semantic.vertex).first!
        let normalSource = geometry.sources(for: SCNGeometrySource.Semantic.normal).first!
        
        let vertexCount = vertexSource.vectorCount
        
        
        let colourList = [vector_float4](repeating: color, count: vertexCount)
        
        let colourData = Data(bytes: colourList, count: colourList.count * MemoryLayout<vector_float4>.size)
        let colourSource = SCNGeometrySource(data: colourData,
                                             semantic: SCNGeometrySource.Semantic.color,
                                             vectorCount: colourList.count,
                                             usesFloatComponents: true,
                                             componentsPerVector: 4,
                                             bytesPerComponent: MemoryLayout<Float>.size,
                                             dataOffset: 0,
                                             dataStride: MemoryLayout<vector_float4>.size)
        let geo = SCNGeometry(sources: [vertexSource,normalSource,colourSource], elements: [indexElement])
        return geo
    }
    
    static func buildPoligon() -> SCNGeometry {
        
        // form
//        let points: [SCNVector3] = [SCNVector3(x: 0, y: 0, z:0),
//                                    SCNVector3(x: -11.8, y: 0.4, z:11.0),
//                                    SCNVector3(x: -15.5, y: 2.1, z:42),
//                                    SCNVector3(x: -16.5, y: 3, z:72.6),
//                                    SCNVector3(x: -6.1, y: 3.4, z:85.1),
//                                    SCNVector3(x: 26.2695923, y: 3.44000006,z : 88.1231613),
//                                    SCNVector3(x : 52.3630257, y : 2.58999991, z : 96.565506),
//                                    SCNVector3(x : 88.8584671, y : 1.30999994, z : 103.776833),
//                                    SCNVector3(x : 125.000931, y : 0.25, z : 112.923729),
//                                    SCNVector3(x : 149.859924, y : -0.170000002, z : 123.300522),
//                                    SCNVector3(x : 172.779846, y : -1.45000005, z : 120.135063),
//                                    SCNVector3(x : 178.773895, y : -2.29999995, z : 102.72168),
//                                    SCNVector3(x : 180.184494, y : -2.50999999, z : 89.3541718),
//                                    SCNVector3(x : 173.661301, y : -2.07999992, z : 77.2176743),
//                                    SCNVector3(x : 87.9763489, y : -1.00999999, z : 43.0943565),
//                                    SCNVector3(x: 0, y: 0, z:0),
//                                    SCNVector3(x: -11.8, y: 0.4, z:11.0)]
        
        // rect
        let points: [SCNVector3] = [SCNVector3(x: 0, y: 0, z:0),
                                    SCNVector3(x: 10, y: 3, z:0),
                                    SCNVector3(x: 10, y: 0, z:10),
                                    SCNVector3(x: 0, y: 2, z:10)]
        
        let radius: Float = 1
        var centerPoints:[SCNVector3] = preparedPoints(points: points, radius: 3*radius)
        centerPoints = preparedPoints(points: centerPoints, radius: radius)
        let geom = buildTube(centerPoints, radius: radius, segmentCount: 4, colour: vector_float4(1.0, 0.0, 0.0, 1.0), secondaryColour: vector_float4(1.0, 1.0, 1.0, 1.0))
        return geom
    }
    
    static func preparedPoints(points: [SCNVector3], radius: Float) -> [SCNVector3] {
        guard points.count > 2 else {
            return points
        }
        var newPoints = [SCNVector3]()
        var prevPoint = points.last!
        for (i, point) in points.enumerated() {
            let nextPoint: SCNVector3
            let j = i + 1
            if j < points.count {
                nextPoint = points[j]
            }
            else {
                nextPoint = points[0]
            }
            
            let beforeL = (point - prevPoint).length()
            let nextL = (nextPoint - point).length()
            let shiftForward = beforeL > 2*radius ? radius : beforeL / 2
            let shiftBackward = nextL > 2*radius ?  radius : nextL / 2
            
            let parallelBefore = (point - prevPoint).normalized() * shiftBackward
            let parallelNext = (nextPoint - point).normalized() * shiftForward
            if parallelBefore.normalized().dot(parallelNext.normalized()) > -0.95 {
                let shift = (parallelNext - parallelBefore) / 4
                newPoints.append(point - parallelBefore)
                newPoints.append(point + shift)
                newPoints.append(point + parallelNext)
            }
            else {
                newPoints.append(point)
            }
            prevPoint = point
        }
        if newPoints.first!.distance(newPoints.last!) > radius / 10 {
            newPoints.append(newPoints.first!)
        }
        return newPoints
    }
    
    static func buildSpiral() -> SCNGeometry {
        
        var centerPoints:[SCNVector3] = []
        
        let angleStep:Float = .pi / 5
        let zStep:Float = 1
        let radius:Float = 1
        
        var angle:Float = 0
        var z:Float = 0
        
        while z < 10 {
            let x = radius * cos(angle)
            let y = radius * sin(angle)
            
            let pt = SCNVector3Make(x, y, z)
            centerPoints.append(pt)
            
            angle += angleStep
            z += zStep
        }
        
        let geom = buildTube(centerPoints, radius: 1, segmentCount: 10, colour: vector_float4(1.0, 0.0, 0.0, 1.0), secondaryColour: vector_float4(1.0, 1.0, 1.0, 1.0))
        return geom
        
    }
    
    
    /**
     Takes the various arrays and makes a geometry object from it
     */
    static func packGeometryTriangles(_ pointsList: [SCNVector3],
                                      indexList: [CInt],
                                      colourList:[vector_float4]?,
                                      normalsList: [SCNVector3]) -> SCNGeometry {
        
        let vertexData = Data(bytes: pointsList, count: pointsList.count * MemoryLayout<vector_float3>.size)
        let vertexSourceNew = SCNGeometrySource(data: vertexData,
                                                semantic: SCNGeometrySource.Semantic.vertex,
                                                vectorCount: pointsList.count,
                                                usesFloatComponents: true,
                                                componentsPerVector: 3,
                                                bytesPerComponent: MemoryLayout<Float>.size,
                                                dataOffset: 0,
                                                dataStride: MemoryLayout<SCNVector3>.size)
        
        let normalData = Data(bytes: normalsList, count: normalsList.count * MemoryLayout<vector_float3>.size)
        let normalSource = SCNGeometrySource(data: normalData,
                                             semantic: SCNGeometrySource.Semantic.normal,
                                             vectorCount: normalsList.count,
                                             usesFloatComponents: true,
                                             componentsPerVector: 3,
                                             bytesPerComponent: MemoryLayout<Float>.size,
                                             dataOffset: 0,
                                             dataStride: MemoryLayout<SCNVector3>.size)
        
        
        let indexData = Data(bytes: indexList, count: indexList.count * MemoryLayout<CInt>.size)
        let indexElement = SCNGeometryElement(
            data: indexData,
            primitiveType: SCNGeometryPrimitiveType.triangles,
            primitiveCount: indexList.count/3,
            bytesPerIndex: MemoryLayout<CInt>.size
        )
        
        if let colourList = colourList {
            let colourData = Data(bytes: colourList, count: colourList.count * MemoryLayout<vector_float4>.size)
            let colourSource = SCNGeometrySource(data: colourData,
                                                 semantic: SCNGeometrySource.Semantic.color,
                                                 vectorCount: colourList.count,
                                                 usesFloatComponents: true,
                                                 componentsPerVector: 4,
                                                 bytesPerComponent: MemoryLayout<Float>.size,
                                                 dataOffset: 0,
                                                 dataStride: MemoryLayout<vector_float4>.size)
            let geo = SCNGeometry(sources: [vertexSourceNew,normalSource,colourSource], elements: [indexElement])
            geo.firstMaterial?.isLitPerPixel = false
            return geo
        } else {
            let geo = SCNGeometry(sources: [vertexSourceNew,normalSource], elements: [indexElement])
            geo.firstMaterial?.isLitPerPixel = false
            return geo
        }
        
    }
    
    
    static func buildTube(_ points: [SCNVector3],
                          radius: Float,
                          segmentCount: Int,
                          colour: vector_float4?,
                          secondaryColour: vector_float4?) -> SCNGeometry {
        
        let isJoined = points[0].distance(points.last!) < radius / 10
        var colour = colour
        if colour == nil {
            colour = vector_float4(0,0,0,1)
        }
        var secondaryColour = secondaryColour
        if secondaryColour == nil {
            secondaryColour = colour
        }
        
        let segmentRotationAngle:Float = 2 * .pi / Float(segmentCount)
        
        var pointsList: [SCNVector3] = []
        var normalsList: [SCNVector3] = []
        var indexList: [CInt] = []
        var colourList:[vector_float4] = []
        
        let lastIndex = isJoined ? points.count - 2 : points.count - 1
        var lastPt:SCNVector3? = points[lastIndex]
        var lastPtsOffset:[SCNVector3]? = nil
        var lastNormalsOffset:[SCNVector3]? = nil
        for (i,pt) in points.enumerated() {
            
            let distanceFromEnd = min(points.count - i, i)
            let radiusMultiplicationFactor: Float = 1 //min(Float(distanceFromEnd) / 4, 1)
            
            if let lastPt = lastPt {
                let along = (pt - lastPt)
                var nextAlong = SCNVector3()
                let j = i + 1
                if j < points.count {
                    let nextPt = points[j]
                    nextAlong = (nextPt - pt)
                }
                else {
                    // if joined figure
                    if isJoined {
                        let nextPt = points[1]
                        let prevPoint = points[0]
                        nextAlong = (nextPt - prevPoint)
                    }
                    else {
                        let nextPt = points[0]
                        nextAlong = (nextPt - pt)
                    }
                }
                
                let parallelBefore = along.normalized()
                let parallelNext = nextAlong.normalized()
                let parallelN = (along + nextAlong).normalized()
                
                
                // find single vector normal to initial vect
                let jointPlane = (parallelBefore - parallelNext).normalized()
                let normalBefore =  (jointPlane - parallelN * jointPlane.dot(parallelBefore)).normalized()
                let coeff = 1 / normalBefore.dot(jointPlane)
                let parallel = parallelN * coeff
                
                let startOffsetVector = SCNVector3(x: 0, y: radius, z: 0)
                let startOffsetVectorGlk = SCNVector3ToGLKVector3(startOffsetVector)
                
                var offsetPoints = [SCNVector3](repeating: SCNVector3Zero, count: segmentCount)
                var offsetNormalVectors = [SCNVector3](repeating: SCNVector3Zero, count: segmentCount)
                
                var rotation:Float = segmentRotationAngle / 2
                for i in 0..<segmentCount {
                    let quart = GLKQuaternionMakeWithAngleAndAxis(rotation, parallel.x, parallel.y, parallel.z)
                    let rotatedOffsetVector = GLKQuaternionRotateVector3(quart, startOffsetVectorGlk)
                    let offsetPoint = pt + SCNVector3FromGLKVector3(rotatedOffsetVector)
                    offsetPoints[i] = offsetPoint
                    offsetNormalVectors[i] = SCNVector3FromGLKVector3(rotatedOffsetVector).normalized()
                    
                    rotation += segmentRotationAngle
                }
                
                if let lastPtsOffset = lastPtsOffset, let lastNormalsOffset = lastNormalsOffset {
                    
                    var prevLastPoint = lastPtsOffset[segmentCount-1]
                    var prevLastNorma = lastNormalsOffset[segmentCount-1]
                    var prevPoint = offsetPoints[segmentCount-1]
                    var prevNorma = offsetNormalVectors[segmentCount-1]
                    for i in 0..<segmentCount {
                        var segmentColour = colour
                        if i == 0 {
                            segmentColour = secondaryColour
                        }
                        
                        let lastPoint = lastPtsOffset[i]
                        let lastNorma = lastNormalsOffset[i]
                        let cPoint = offsetPoints[i]
                        let cNorma = offsetNormalVectors[i]
                        
                        indexList.append(CInt(pointsList.count))
                        pointsList.append(prevPoint)
                        colourList.append(segmentColour!)
                        normalsList.append(prevNorma * 2)
                        
                        indexList.append(CInt(pointsList.count))
                        pointsList.append(prevLastPoint)
                        colourList.append(segmentColour!)
                        normalsList.append(prevLastNorma * 2)
                        
                        indexList.append(CInt(pointsList.count))
                        pointsList.append(lastPoint)
                        colourList.append(segmentColour!)
                        normalsList.append(lastNorma * 2)
                        
                        
                        indexList.append(CInt(pointsList.count))
                        pointsList.append(lastPoint)
                        colourList.append(segmentColour!)
                        normalsList.append(lastNorma)
                        
                        indexList.append(CInt(pointsList.count))
                        pointsList.append(cPoint)
                        colourList.append(segmentColour!)
                        normalsList.append(cNorma)
                        
                        indexList.append(CInt(pointsList.count))
                        pointsList.append(prevPoint)
                        colourList.append(segmentColour!)
                        normalsList.append(prevNorma)
                        
                        
                        prevLastPoint = lastPoint
                        prevLastNorma = lastNorma
                        prevPoint = cPoint
                        prevNorma = cNorma
                    }
                }
                lastPtsOffset = offsetPoints
                lastNormalsOffset = offsetNormalVectors
            }
            lastPt = pt
        }
        
        return GeometryBuildTools.packGeometryTriangles(pointsList, indexList: indexList, colourList: colourList, normalsList: normalsList)
    }
    
}

let scene = SCNScene()
// create the text geometry
let hello = SCNText(string: "Hello\nworld!", extrusionDepth: 9)
hello.chamferRadius = 0.5
hello.flatness = 0.05
hello.font = UIFont(name: "Superclarendon-Black", size: 12)
let helloNode = SCNNode(geometry: hello)

/*: 
 ## Positioning the text geometry
 
 SCNText geometries have their origin in the lower-left corner.
 We want to move the pivot point of the text geometry to its centre.
 
 First, we get the bounding box of the geometry, expressed as two SCNVector3 vectors, describing the distance from the pivot point to the lower-left-back corner and upper-right-front corner respectively
 */
var scnMin = SCNVector3Zero // this will hold distance from pivot to lower-left-back...
var scnMax = SCNVector3Zero // distance from pivot to upper-right-front
//helloNode.__getBoundingBoxMin(&scnMin, max: &scnMax) // pass min and max as in-out parameters. in beta 4 getBoundingBox now has 2 underscores in its method name?

/*:
 SCNVector3 vectors do not come with `+ - * /` operators as standard. A common approach in Swift is to write SCNvector3 operators. Instead, here we'll bridge to simd, which does handle `+ - *` operators (though not scalar `/` ). simd is imported by SceneKit and SpriteKit.
 
 We need to move the pivot by half the total length - the current position of the pivot.
 
 ```
 halfTotalLength = (min + max) * 0.5 
 // nb scalar / not supported
 pivotPosition = min
 
 translation = (( min + max) * 0.5 ) - min 
 == (max - min) * 0.5
 ```

 */
let minScnFloat = float3.init(scnMin)
let maxScnFloat = float3.init(scnMax)
let translation = (maxScnFloat - minScnFloat) * 0.5
helloNode.pivot = SCNMatrix4MakeTranslation(translation.x, translation.y, translation.z)
 
//helloNode.scale = SCNVector3(0.2, 0.2, 0.2)
//scene.rootNode.addChildNode(helloNode)

// The cube will be an underscore for the "Hello"
let cube = SCNNode(geometry: SCNBox(width: 30, height: 2, length: 8, chamferRadius: 0.5))
cube.position = SCNVector3(-4, 3, 0)
//scene.rootNode.addChildNode(cube)





///Flat Bezier Path extruded

let path = UIBezierPath()
path.move(to: CGPoint.zero)
path.addQuadCurve(to: CGPoint(x: 100, y: 0), controlPoint: CGPoint(x: 50, y: 200))
path.addLine(to: CGPoint(x: 98, y: 0))
path.addQuadCurve(to: CGPoint(x: 2, y: 0), controlPoint: CGPoint(x: 50, y: 198))

// Tweak for a smoother shape (lower is smoother)
path.flatness = 0.25

// Make a 3D extruded shape from the path
let shape = SCNShape(path: path, extrusionDepth: 10)
shape.firstMaterial?.diffuse.contents = UIColor(white: 1, alpha: 1)

// And place it in the scene
let shapeNode = SCNNode(geometry: shape)
shapeNode.pivot = SCNMatrix4MakeTranslation(50, 0, 0)
shapeNode.eulerAngles.y = Float(-M_PI_4)
//scene.rootNode.addChildNode(shapeNode)


let pathFill = UIBezierPath()
path.move(to: CGPoint.zero)
path.addQuadCurve(to: CGPoint(x: 99, y: 0), controlPoint: CGPoint(x: 50, y: 199))

// Tweak for a smoother shape (lower is smoother)
path.flatness = 0.25

// Make a 3D extruded shape from the path
let shapeFill = SCNShape(path: path, extrusionDepth: 1)
shapeFill.firstMaterial?.diffuse.contents = UIColor(red: 0, green: 1, blue: 0, alpha: 0.2)

// And place it in the scene
let shapeNodeFill = SCNNode(geometry: shapeFill)
shapeNodeFill.pivot = SCNMatrix4MakeTranslation(50, 0, 0)
shapeNodeFill.eulerAngles.y = Float(-M_PI_4)
//scene.rootNode.addChildNode(shapeNodeFill)

let polyline = [CGPoint(x:0, y:0), CGPoint(x:0, y:100)]
//
let spiral = GeometryBuildTools.buildPoligon() //SCNSphere(radius: 0.5)
//spiral.firstMaterial?.diffuse.contents = UIColor.red
//spiral.firstMaterial?.selfIllumination.contents = UIColor.red

let parent = SCNNode()
let sphereNode = SCNNode(geometry: spiral)
sphereNode.opacity = 0.4
parent.addChildNode(sphereNode)
scene.rootNode.addChildNode(parent)




/*: 
 ## Add materials
 
 Set up two materials with physics-based rendering (PBR). Note that the PBR properties won't show up when running the playground in Xcode (I think because PBR requires Metal, and the inline iOS Simulator is backed by OpenGL, not Metal). If you want to experiment with PBR on the Mac, you need to have MacOS Sierra installed, and set up a new Playground as a MacOS Playground and port this code over.
 */
let dullMat = SCNMaterial()
dullMat.diffuse.contents = #colorLiteral(red: 0.7602152824, green: 0.7601925135, blue: 0.7602053881, alpha: 1)
let metalMat = SCNMaterial()
metalMat.diffuse.contents = #colorLiteral(red: 1.0, green: 0.498039215803146, blue: 0.756862759590149, alpha: 1.0)

dullMat.lightingModel = SCNMaterial.LightingModel.physicallyBased
dullMat.roughness.contents = #colorLiteral(red: 1, green: 1, blue: 1, alpha: 1)
dullMat.metalness.contents = #colorLiteral(red: 0.1956433058, green: 0.2113749981, blue: 0.2356699705, alpha: 1)

metalMat.lightingModel = SCNMaterial.LightingModel.physicallyBased
metalMat.roughness.contents = #colorLiteral(red: 0.5296475887, green: 0.5296317339, blue: 0.5296407342, alpha: 1)
//add some noise to the metal texture using SpriteKit. You could really go to town with the new GameplayKit noise features
metalMat.metalness.contents = SKTexture(noiseWithSmoothness: 0.8, size: CGSize(width: 500, height: 500), grayscale: true).cgImage()

cube.geometry?.materials = [dullMat]
// SCNText geometries have up to 5 materials: front, back, sides, front champfer, back champfer
// Apply the metallic texture to only the front face and front champfers:
hello.materials = [metalMat, dullMat, dullMat, metalMat, dullMat]

let light = SCNLight()
light.type = SCNLight.LightType.omni
let lightNode = SCNNode()
lightNode.light = light
lightNode.position = SCNVector3(8,12,15)
//scene.rootNode.addChildNode(lightNode)

/*:
 ## Set up the live view
 
 In Xcode on the Mac, you have to supply a size for the view. On iPad, this isn't necessary
 */

//let view = SCNView() //iPad version
let view = SCNView(frame: CGRect(x: 0, y: 0, width: 600, height: 600)) //Xcode version
view.allowsCameraControl = true
view.autoenablesDefaultLighting = true
view.showsStatistics = true
view.scene = scene
view.backgroundColor = #colorLiteral(red: 0.0588235296308994, green: 0.180392161011696, blue: 0.24705882370472, alpha: 1.0)
PlaygroundPage.current.liveView = view
