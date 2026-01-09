import SwiftUI

private enum Constants {
    // Information Theory
    static let bitsPerMolecule: UInt8 = 8
    static let energyPerBit: Double = 4
    static let photonEnergy: Double = 15
    static let complementarityBonus: Double = 3

    // Thermodynamics
    static let activationTemperature: Double = 1.5
    static let entropyRate: Double = 0.08
    static let maintenanceCost: Double = 0.02
    static let entropyDeathThreshold: Double = 150
    static let entropyRepairFactor: Double = 0.94
    static let entropyRepairCost: Double = 0.8
    static let entropyChildFactor: Double = 0.8
    static let entropyParentFactor: Double = 1.2

    // Chemistry
    static let reactionRate: Int = 8
    static let catalysisThreshold: Int = 3
    static let minRepairMolecules: Int = 5
    static let minRepairEnergy: Double = 5

    // Physics
    static let lightAttenuation: Double = 3.0
    static let pressureGradient: Double = 8
    static let temperatureGradient: Double = 12
    static let minTemperature: Double = 0
    static let diffusionRate: Double = 0.008
    static let gravity: Double = 0.03
    static let friction: Double = 0.96

    // Replication
    static let cloneEnergyThreshold: Double = 350
    static let mutationRate: Double = 0.02
    static let duplicationRate: Double = 0.01
    static let deletionRate: Double = 0.01

    // Initial Conditions
    static let initialEnergy: Double = 120
    static let initialEntropy: Double = 10

    // Grid
    static let debrisGridSize: Double = 18
    static let nutrientGridSize: Double = 25
}

final class SimulationModel: ObservableObject {
    @Published var population: Int = 0
    @Published var averageEnergy: Int = 0
    @Published var averageEntropy: Int = 0
    @Published var speed: Int = 1

    private var environment = Environment(width: 1, height: 1)
    private var cells: [Cell] = []
    private var frameCounter = 0

    func updateSize(_ size: CGSize) {
        guard size.width > 0, size.height > 0 else { return }
        environment = Environment(width: size.width, height: size.height)
    }

    func spawn(count: Int) {
        guard environment.width > 0, environment.height > 0 else { return }
        for _ in 0..<count {
            let genomeLength = 6 + Int.random(in: 0..<8)
            let genome = (0..<genomeLength).map { _ in UInt8.random(in: 0...255) }
            let cell = Cell(genome: genome,
                            x: Double.random(in: 0..<environment.width),
                            y: Double.random(in: 0..<environment.height))
            cells.append(cell)
        }
        updateStats()
    }

    func reset() {
        cells.removeAll()
        environment = Environment(width: environment.width, height: environment.height)
        updateStats()
    }

    func stepFrame() {
        guard environment.width > 0, environment.height > 0 else { return }
        for _ in 0..<max(speed, 1) {
            let currentTime = Date().timeIntervalSince1970
            environment.step()
            resolveCollisions()

            var babies: [Cell] = []
            for cell in cells {
                cell.timestep(environment: environment, time: currentTime)
                if cell.alive && cell.energy > Constants.cloneEnergyThreshold {
                    if cells.count + babies.count < 500 {
                        babies.append(cell.clone())
                    }
                }
            }
            cells = cells.filter { $0.alive }
            cells.append(contentsOf: babies)
        }

        frameCounter += 1
        if frameCounter % 6 == 0 {
            updateStats()
        }
    }

    func draw(in context: inout GraphicsContext, size: CGSize) {
        guard size.width > 0, size.height > 0 else { return }
        if environment.width != size.width || environment.height != size.height {
            updateSize(size)
        }

        let rect = CGRect(origin: .zero, size: size)
        let gradient = Gradient(colors: [
            Color(red: 0.06, green: 0.65, blue: 0.91),
            Color(red: 0.12, green: 0.23, blue: 0.54),
            Color(red: 0.06, green: 0.09, blue: 0.16),
            Color(red: 0.1, green: 0.08, blue: 0.08)
        ])

        context.fill(Path(rect), with: .linearGradient(gradient,
                                                      startPoint: CGPoint(x: 0, y: 0),
                                                      endPoint: CGPoint(x: 0, y: size.height)))

        drawDebris(in: &context)

        for cell in cells {
            let energyColor = energyColor(for: cell.energy)
            let entropyRatio = min(1, cell.entropy / Constants.entropyDeathThreshold)
            context.opacity = 1.0 - (entropyRatio * 0.4)

            let center = CGPoint(x: cell.position.x, y: cell.position.y)
            let bodyRect = CGRect(x: center.x - cell.radius,
                                  y: center.y - cell.radius,
                                  width: cell.radius * 2,
                                  height: cell.radius * 2)
            context.fill(Path(ellipseIn: bodyRect), with: .color(energyColor))

            let membraneCount = cell.membrane.count
            if membraneCount > 0 {
                var strokeStyle = StrokeStyle(lineWidth: 1 + Double(membraneCount) / 8)
                context.stroke(Path(ellipseIn: bodyRect), with: .color(.white.opacity(0.3)), style: strokeStyle)
            }

            context.opacity = 1.0
            let nucleusRect = CGRect(x: center.x - 3, y: center.y - 3, width: 6, height: 6)
            context.fill(Path(ellipseIn: nucleusRect), with: .color(cell.strategyColor()))
            context.stroke(Path(ellipseIn: nucleusRect), with: .color(.black.opacity(0.5)), style: StrokeStyle(lineWidth: 1))

            let ringRect = CGRect(x: center.x - (cell.radius + 2),
                                  y: center.y - (cell.radius + 2),
                                  width: (cell.radius + 2) * 2,
                                  height: (cell.radius + 2) * 2)
            context.stroke(Path(ellipseIn: ringRect), with: .color(cell.genomeColor), style: StrokeStyle(lineWidth: 0.5))
        }
        context.opacity = 1.0
    }

    private func energyColor(for energy: Double) -> Color {
        if energy < 100 {
            let p = max(0, energy / 100)
            return Color(red: 1.0, green: p, blue: 0)
        }
        let p = min(1, (energy - 100) / 200)
        return Color(red: 1 - p, green: 1.0, blue: 0)
    }

    private func drawDebris(in context: inout GraphicsContext) {
        for (key, molecules) in environment.debris.grid {
            let total = molecules.values.reduce(0, +)
            let size = min(Constants.debrisGridSize - 2, max(2, Double(total) / 3))
            let rect = CGRect(x: Double(key.x) * Constants.debrisGridSize + Constants.debrisGridSize / 2 - size / 2,
                              y: Double(key.y) * Constants.debrisGridSize + Constants.debrisGridSize / 2 - size / 2,
                              width: size,
                              height: size)
            context.fill(Path(rect), with: .color(Color(red: 0.63, green: 0.55, blue: 0.7, opacity: 0.4)))
        }
    }

    private func resolveCollisions() {
        guard cells.count > 1 else { return }
        for i in 0..<(cells.count - 1) {
            for j in (i + 1)..<cells.count {
                let c1 = cells[i]
                let c2 = cells[j]
                let dx = c1.position.x - c2.position.x
                let dy = c1.position.y - c2.position.y
                let distSq = dx * dx + dy * dy
                let minDist = c1.radius + c2.radius

                if distSq > 0 && distSq < minDist * minDist {
                    let dist = sqrt(distSq)
                    let overlap = minDist - dist
                    let nx = dx / dist
                    let ny = dy / dist
                    let force = overlap * 0.05

                    c1.velocity.dx += nx * force
                    c1.velocity.dy += ny * force
                    c2.velocity.dx -= nx * force
                    c2.velocity.dy -= ny * force
                }
            }
        }
    }

    private func updateStats() {
        population = cells.count
        guard !cells.isEmpty else {
            averageEnergy = 0
            averageEntropy = 0
            return
        }
        let totalEnergy = cells.reduce(0) { $0 + $1.energy }
        let totalEntropy = cells.reduce(0) { $0 + $1.entropy }
        averageEnergy = Int(totalEnergy / Double(cells.count))
        averageEntropy = Int(totalEntropy / Double(cells.count))
    }
}

private struct MoleculeProperties {
    let complexity: Int
    let stability: Int
    let reactivity: Int
    let polarity: Int
    let conjugation: Int
    let size: Int
}

private func popcount(_ value: UInt8) -> Int {
    var count = 0
    var n = value
    while n > 0 {
        if n & 1 == 1 { count += 1 }
        n >>= 1
    }
    return count
}

private func countAlternatingBits(_ value: UInt8) -> Int {
    var alternating = 0
    for i in 0..<7 {
        let bit1 = (value >> i) & 1
        let bit2 = (value >> (i + 1)) & 1
        if bit1 != bit2 {
            alternating += 1
        }
    }
    return alternating
}

private func getMoleculeProperties(_ molecule: UInt8) -> MoleculeProperties {
    let bits = popcount(molecule)
    return MoleculeProperties(
        complexity: bits,
        stability: 8 - abs(bits - 4),
        reactivity: popcount(molecule & 0x0F),
        polarity: popcount(molecule & 0xAA) - popcount(molecule & 0x55),
        conjugation: countAlternatingBits(molecule),
        size: bits / 2
    )
}

private struct ReactionResult {
    let products: [UInt8]
    let energy: Double
}

private func attemptReaction(_ mA: UInt8, _ mB: UInt8) -> ReactionResult {
    let pathways: [[UInt8]] = [
        [mA ^ mB, mA & mB],
        [mA | mB, mA & ~mB],
        [UInt8((UInt16(mA) + UInt16(mB)) >> 1), ((mA & 1) == 1 ? 0xFF : 0x00)],
        [~(mA ^ mB), mA ^ 0xFF]
    ]

    var best = ReactionResult(products: [], energy: -Double.infinity)
    for products in pathways {
        let reactantBits = popcount(mA) + popcount(mB)
        let productBits = products.reduce(0) { $0 + popcount($1) }
        let orderEnergy = Double(reactantBits - productBits) * Constants.energyPerBit
        let complementarity = products.count >= 2 ? popcount(products[0] & products[1]) : 0
        let totalEnergy = orderEnergy + Double(complementarity) * Constants.complementarityBonus
        let filtered = products.filter { $0 > 0 }
        if totalEnergy > best.energy {
            best = ReactionResult(products: filtered, energy: totalEnergy)
        }
    }
    return best
}

private struct PhysicsSample {
    let lightIntensity: Double
    let pressure: Double
    let temperature: Double
    let turbulence: Double
}

private final class Environment {
    let width: Double
    let height: Double
    var nutrientField: [Double]
    let nutrientGridWidth: Int
    let nutrientGridHeight: Int
    let debris: DebrisManager

    init(width: Double, height: Double) {
        self.width = width
        self.height = height
        nutrientGridWidth = Int(ceil(width / Constants.nutrientGridSize))
        nutrientGridHeight = Int(ceil(height / Constants.nutrientGridSize))
        nutrientField = Array(repeating: 1.0, count: nutrientGridWidth * nutrientGridHeight)
        debris = DebrisManager(width: width, height: height)
    }

    func getPhysics(x: Double, y: Double, time: Double) -> PhysicsSample {
        let depth = y / height
        let lightIntensity = exp(-depth * Constants.lightAttenuation)
        let pressure = depth * Constants.pressureGradient
        let temperature = max(Constants.minTemperature, 20 - depth * Constants.temperatureGradient)
        let turbulence = max(0, (1 - depth) * sin(x / 50 + time) * 0.1)
        return PhysicsSample(lightIntensity: lightIntensity,
                             pressure: pressure,
                             temperature: temperature,
                             turbulence: turbulence)
    }

    func tryGetNutrient(x: Double, y: Double, traits: CellTraits, time: Double) -> UInt8? {
        let physics = getPhysics(x: x, y: y, time: time)
        let depth = y / height

        let gx = min(Int(x / Constants.nutrientGridSize), nutrientGridWidth - 1)
        let gy = min(Int(y / Constants.nutrientGridSize), nutrientGridHeight - 1)
        let idx = gy * nutrientGridWidth + gx

        if nutrientField[idx] < 0.1 { return nil }

        if Double.random(in: 0..<1) < physics.lightIntensity * 0.15 * traits.uptakeRate {
            nutrientField[idx] -= 0.01
            return 0x03
        }

        if physics.pressure > 5 && Double.random(in: 0..<1) < 0.08 * traits.uptakeRate {
            nutrientField[idx] -= 0.01
            return 0xFE
        }

        let debrisDensity = debris.getDensity(x: x, y: y)
        if Double.random(in: 0..<1) < debrisDensity * 0.02 * traits.scavengeAbility {
            return debris.scavenge(x: x, y: y)
        }

        if depth > 0.6 && Double.random(in: 0..<1) < 0.05 * traits.uptakeRate {
            nutrientField[idx] -= 0.015
            return 0xC6
        }

        return nil
    }

    func step() {
        for i in nutrientField.indices {
            let gy = i / nutrientGridWidth
            let depth = Double(gy) / Double(nutrientGridHeight)
            let regenRate = exp(-depth * 2) * Constants.diffusionRate
            nutrientField[i] = min(1.0, nutrientField[i] + regenRate)
        }
        debris.step()
    }
}

private final class DebrisManager {
    let width: Double
    let height: Double
    struct GridKey: Hashable {
        let x: Int
        let y: Int
    }
    var grid: [GridKey: [UInt8: Int]] = [:]

    init(width: Double, height: Double) {
        self.width = width
        self.height = height
    }

    private func key(x: Double, y: Double) -> GridKey {
        let gx = Int(x / Constants.debrisGridSize)
        let gy = Int(y / Constants.debrisGridSize)
        return GridKey(x: gx, y: gy)
    }

    func add(x: Double, y: Double, molecules: [UInt8: Int]) {
        let k = key(x: x, y: y)
        if grid[k] == nil {
            grid[k] = [:]
        }
        for (mol, count) in molecules {
            grid[k, default: [:]][mol, default: 0] += count
        }
    }

    func scavenge(x: Double, y: Double) -> UInt8? {
        let k = key(x: x, y: y)
        guard var cellDebris = grid[k], !cellDebris.isEmpty else { return nil }
        let types = Array(cellDebris.keys)
        guard let type = types.randomElement() else { return nil }
        if let amount = cellDebris[type], amount > 0 {
            cellDebris[type] = amount - 1
            if cellDebris[type] == 0 {
                cellDebris[type] = nil
            }
            if cellDebris.isEmpty {
                grid[k] = nil
            } else {
                grid[k] = cellDebris
            }
            return type
        }
        return nil
    }

    func getDensity(x: Double, y: Double) -> Double {
        let k = key(x: x, y: y)
        guard let cellDebris = grid[k] else { return 0 }
        let total = cellDebris.values.reduce(0, +)
        return min(1.0, Double(total) / 50)
    }

    func step() {
        let maxGy = Int(height / Constants.debrisGridSize) - 1
        let sortedKeys = grid.keys.sorted { lhs, rhs in
            lhs.y > rhs.y
        }

        for key in sortedKeys {
            if key.y >= maxGy { continue }

            if Double.random(in: 0..<1) < Constants.gravity {
                let targetKey = GridKey(x: key.x, y: key.y + 1)
                if grid[targetKey] == nil {
                    grid[targetKey] = [:]
                }
                let source = grid[key] ?? [:]
                for (mol, count) in source {
                    grid[targetKey, default: [:]][mol, default: 0] += count
                }
                grid[key] = nil
            }
        }

        if Double.random(in: 0..<1) < Constants.entropyRate * 0.5 {
            guard let randomKey = grid.keys.randomElement() else { return }
            guard var cellDebris = grid[randomKey] else { return }
            guard let randomType = cellDebris.keys.randomElement() else { return }
            cellDebris[randomType, default: 0] -= 1
            if cellDebris[randomType] ?? 0 <= 0 {
                cellDebris[randomType] = nil
            }
            grid[randomKey] = cellDebris.isEmpty ? nil : cellDebris
        }
    }
}

private struct CellTraits {
    let uptakeRate: Double
    let permeability: Double
    let scavengeAbility: Double
    let metabolicRate: Double
}

private struct CellStats {
    var lightEnergy: Double = 0
    var catalyticEnergy: Double = 0
    var scavengeEnergy: Double = 0
}

private final class Cell {
    let genome: [UInt8]
    private var ribosomePos = 0
    private(set) var cytoplasm: [UInt8: Int] = [:]
    private(set) var membrane: [UInt8: Int] = [:]
    var energy: Double = Constants.initialEnergy
    var entropy: Double = Constants.initialEntropy
    var position: CGPoint
    var velocity: CGVector
    var age = 0
    var alive = true
    var stats = CellStats()
    let radius: Double
    let genomeColor: Color

    private var cachedTraits: CellTraits?
    private var cachedTraitsAge: Int = 0

    init(genome: [UInt8], x: Double, y: Double) {
        self.genome = genome
        position = CGPoint(x: x, y: y)
        velocity = CGVector(dx: Double.random(in: -1...1), dy: Double.random(in: -1...1))
        radius = 4 + Double(genome.count) / 6
        genomeColor = Color(hue: Double(genome.first ?? 0) / 255.0, saturation: 0.6, brightness: 0.5)
    }

    var traits: CellTraits {
        if let cachedTraits, age % 20 != 0 { return cachedTraits }
        let avgBits = Double(genome.reduce(0) { $0 + popcount($1) }) / Double(genome.count)
        let diversity = Double(Set(genome).count) / Double(genome.count)
        let symmetry = Double(genome.enumerated().filter { index, gene in
            gene == genome[genome.count - 1 - index]
        }.count) / Double(genome.count)

        let traits = CellTraits(
            uptakeRate: min(1, avgBits / 6),
            permeability: diversity * 0.8,
            scavengeAbility: (1 - symmetry) * 0.7,
            metabolicRate: 1 / log(Double(genome.count) + 2)
        )
        cachedTraits = traits
        cachedTraitsAge = age
        return traits
    }

    func timestep(environment: Environment, time: Double) {
        guard alive else { return }
        age += 1
        let physics = environment.getPhysics(x: position.x, y: position.y, time: time)

        entropy += Constants.entropyRate
        let maintenanceCost = entropy * Constants.maintenanceCost
        energy -= maintenanceCost

        expressionStep()

        if canRepair() {
            entropy *= Constants.entropyRepairFactor
            energy -= Constants.entropyRepairCost
        }

        photoreceptionStep(intensity: physics.lightIntensity)

        if let nutrient = environment.tryGetNutrient(x: position.x, y: position.y, traits: traits, time: time) {
            addToPool(&cytoplasm, molecule: nutrient, amount: 1)
        }

        for _ in 0..<Constants.reactionRate {
            catalyticReactionStep()
        }

        if age % 5 == 0 {
            motorStep(environment: environment, time: time)
        }

        velocity.dx += (Double.random(in: 0..<1) - 0.5) * 0.1 + physics.turbulence
        velocity.dy += (Double.random(in: 0..<1) - 0.5) * 0.1 + Constants.gravity * 0.3

        position.x += velocity.dx
        position.y += velocity.dy
        velocity.dx *= Constants.friction
        velocity.dy *= Constants.friction

        if position.x < radius {
            position.x = radius
            velocity.dx *= -0.8
        }
        if position.x > environment.width - radius {
            position.x = environment.width - radius
            velocity.dx *= -0.8
        }
        if position.y < radius {
            position.y = radius
            velocity.dy *= -0.8
        }
        if position.y > environment.height - radius {
            position.y = environment.height - radius
            velocity.dy *= -0.8
        }

        if entropy > Constants.entropyDeathThreshold || energy < 0 {
            die(environment: environment)
        }
    }

    private func expressionStep() {
        guard !genome.isEmpty else { return }
        let gene = genome[ribosomePos]
        let totalMolecules = cytoplasm.values.reduce(0, +)

        if totalMolecules < 80 && energy >= 1 {
            addToPool(&cytoplasm, molecule: gene, amount: 1)
            energy -= 1

            let props = getMoleculeProperties(gene)
            if Double.random(in: 0..<1) < traits.permeability && props.polarity < 2 {
                addToPool(&cytoplasm, molecule: gene, amount: -1)
                addToPool(&membrane, molecule: gene, amount: 1)
            }
        }

        ribosomePos = (ribosomePos + 1) % genome.count
    }

    private func canRepair() -> Bool {
        let mols = poolList(cytoplasm)
        return mols.count > Constants.minRepairMolecules && energy > Constants.minRepairEnergy
    }

    private func catalyticReactionStep() {
        let mols = poolList(cytoplasm)
        guard mols.count >= 3 else { return }
        let catalystIndex = Int.random(in: 0..<mols.count)
        let catalyst = mols[catalystIndex]
        let availableIndices = mols.indices.filter { $0 != catalystIndex }
        guard availableIndices.count >= 2 else { return }
        let sub1Index = availableIndices.randomElement() ?? catalystIndex
        let sub2Indices = availableIndices.filter { $0 != sub1Index }
        guard let sub2Index = sub2Indices.randomElement() else { return }
        let sub1 = mols[sub1Index]
        let sub2 = mols[sub2Index]

        let catalyticFit = popcount(catalyst & sub1 & sub2)
        if catalyticFit >= Constants.catalysisThreshold {
            let reaction = attemptReaction(sub1, sub2)
            if energy >= max(0, -reaction.energy * 0.1) {
                addToPool(&cytoplasm, molecule: sub1, amount: -1)
                addToPool(&cytoplasm, molecule: sub2, amount: -1)
                for product in reaction.products {
                    addToPool(&cytoplasm, molecule: product, amount: 1)
                }
                energy += reaction.energy
                if reaction.energy > 0 {
                    stats.catalyticEnergy += reaction.energy
                }
            }
        }
    }

    private func photoreceptionStep(intensity: Double) {
        if Double.random(in: 0..<1) > intensity || Double.random(in: 0..<1) > 0.5 {
            return
        }
        let allMolecules = poolList(membrane) + poolList(cytoplasm)
        guard let target = allMolecules.randomElement() else { return }
        let props = getMoleculeProperties(target)
        let absorptionChance = Double(props.conjugation) / 7.0
        if Double.random(in: 0..<1) < absorptionChance {
            let flip = UInt8(1 << Int.random(in: 0..<Int(Constants.bitsPerMolecule)))
            let excited = (target ^ flip) & 0xFF
            let energyGain = Constants.photonEnergy * (Double(props.conjugation) / 7.0)

            if membrane[target, default: 0] > 0 {
                addToPool(&membrane, molecule: target, amount: -1)
                addToPool(&membrane, molecule: excited, amount: 1)
            } else {
                addToPool(&cytoplasm, molecule: target, amount: -1)
                addToPool(&cytoplasm, molecule: excited, amount: 1)
            }

            energy += energyGain
            stats.lightEnergy += energyGain
        }
    }

    private func motorStep(environment: Environment, time: Double) {
        let mems = poolList(membrane)
        let fuels = poolList(cytoplasm)
        guard !mems.isEmpty, !fuels.isEmpty, energy >= 3 else { return }

        let gradients = sampleGradients(environment: environment, time: time)
        var netForce = CGVector(dx: 0, dy: 0)
        let sampleSize = min(5, mems.count)

        for _ in 0..<sampleSize {
            guard let mem = mems.randomElement() else { continue }
            for gradient in gradients {
                let affinity = Double(popcount(mem & gradient.moleculeType)) / 8.0
                if affinity > 0.4 {
                    let memPolarity = getMoleculeProperties(mem).polarity
                    let gradPolarity = getMoleculeProperties(gradient.moleculeType).polarity
                    let attraction: Double = (memPolarity.signum() == gradPolarity.signum()) ? 1 : -1
                    netForce.dx += gradient.dx * affinity * attraction * 0.5
                    netForce.dy += gradient.dy * affinity * attraction * 0.5
                }
            }
        }

        if abs(netForce.dx) > 0.1 || abs(netForce.dy) > 0.1 {
            let fuel = fuels[0]
            let fuelPower = Double(popcount(fuel)) * 0.15
            velocity.dx += netForce.dx * fuelPower
            velocity.dy += netForce.dy * fuelPower
            addToPool(&cytoplasm, molecule: fuel, amount: -1)
            energy -= 3
        }
    }

    private func sampleGradients(environment: Environment, time: Double) -> [GradientSample] {
        let sampleDist: Double = 30
        let directions: [CGVector] = [
            CGVector(dx: sampleDist, dy: 0),
            CGVector(dx: -sampleDist, dy: 0),
            CGVector(dx: 0, dy: sampleDist),
            CGVector(dx: 0, dy: -sampleDist)
        ]
        var gradients: [GradientSample] = []

        for dir in directions {
            let sampleX = max(0, min(environment.width - 1, position.x + dir.dx))
            let sampleY = max(0, min(environment.height - 1, position.y + dir.dy))
            let physics = environment.getPhysics(x: sampleX, y: sampleY, time: time)
            let density = environment.debris.getDensity(x: sampleX, y: sampleY)

            if physics.lightIntensity > 0.1 {
                gradients.append(GradientSample(moleculeType: 0x03,
                                               concentration: physics.lightIntensity,
                                               dx: dir.dx.signum(),
                                               dy: dir.dy.signum()))
            }

            if density > 0.1 {
                gradients.append(GradientSample(moleculeType: 0xC6,
                                               concentration: density,
                                               dx: dir.dx.signum(),
                                               dy: dir.dy.signum()))
            }
        }

        return gradients
    }

    func clone() -> Cell {
        var newGenome = genome
        for i in newGenome.indices {
            if Double.random(in: 0..<1) < Constants.mutationRate {
                let bit = UInt8(1 << Int.random(in: 0..<Int(Constants.bitsPerMolecule)))
                newGenome[i] ^= bit
            }
        }

        if Double.random(in: 0..<1) < Constants.duplicationRate && newGenome.count < 20 {
            let index = Int.random(in: 0..<newGenome.count)
            newGenome.insert(newGenome[index], at: index)
        }

        if Double.random(in: 0..<1) < Constants.deletionRate && newGenome.count > 3 {
            newGenome.remove(at: Int.random(in: 0..<newGenome.count))
        }

        let child = Cell(genome: newGenome, x: position.x, y: position.y)
        let splitEnergy = energy / 2
        child.energy = splitEnergy
        energy = splitEnergy
        child.entropy = entropy * Constants.entropyChildFactor
        entropy *= Constants.entropyParentFactor
        return child
    }

    private func die(environment: Environment) {
        alive = false
        var corpse = cytoplasm
        for (mol, count) in membrane {
            corpse[mol, default: 0] += count
        }
        environment.debris.add(x: position.x, y: position.y, molecules: corpse)
    }

    func strategyColor() -> Color {
        let total = stats.lightEnergy + stats.catalyticEnergy + stats.scavengeEnergy
        if total < 5 { return Color.gray }
        let light = stats.lightEnergy / total
        let catalytic = stats.catalyticEnergy / total
        let scavenge = stats.scavengeEnergy / total

        if light > catalytic && light > scavenge { return Color.yellow }
        if catalytic > light && catalytic > scavenge { return Color.green }
        if scavenge > light && scavenge > catalytic { return Color.purple }
        return Color.white
    }

    private func addToPool(_ pool: inout [UInt8: Int], molecule: UInt8, amount: Int) {
        pool[molecule, default: 0] += amount
        if pool[molecule, default: 0] <= 0 {
            pool[molecule] = nil
        }
    }

    private func poolList(_ pool: [UInt8: Int]) -> [UInt8] {
        var list: [UInt8] = []
        for (mol, count) in pool {
            if count > 0 {
                list.append(contentsOf: Array(repeating: mol, count: count))
            }
        }
        return list
    }
}

private struct GradientSample {
    let moleculeType: UInt8
    let concentration: Double
    let dx: Double
    let dy: Double
}
