//
//  ContentView.swift
//  test2
//
//  Created by Thomas Holeva on 1/9/26.
//

import SwiftUI

struct ContentView: View {
    @StateObject private var simulation = SimulationModel()
    private let timer = Timer.publish(every: 1.0 / 60.0, on: .main, in: .common).autoconnect()

    var body: some View {
        GeometryReader { proxy in
            ZStack(alignment: .top) {
                Canvas { context, size in
                    simulation.draw(in: &context, size: size)
                }
                .onAppear {
                    simulation.updateSize(proxy.size)
                    if simulation.population == 0 {
                        simulation.spawn(count: 10)
                    }
                }
                .onChange(of: proxy.size) { _, newSize in
                    simulation.updateSize(newSize)
                }
                .ignoresSafeArea()

                VStack(spacing: 8) {
                    HStack(alignment: .center, spacing: 12) {
                        StatView(label: "POP", value: "\(simulation.population)", color: .emerald)
                        StatView(label: "AVG NRG", value: "\(simulation.averageEnergy)", color: .amber)
                        StatView(label: "AVG ENT", value: "\(simulation.averageEntropy)", color: .red)
                        Divider()
                            .frame(height: 16)
                        HStack(spacing: 8) {
                            LegendItem(color: .yellow, label: "High-Energy")
                            LegendItem(color: .green, label: "Catalytic")
                            LegendItem(color: .purple, label: "Scavenger")
                        }
                    }
                    .padding(.horizontal, 12)
                    .padding(.vertical, 8)
                    .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 10))

                    HStack(spacing: 8) {
                        Button("+ Random") {
                            simulation.spawn(count: 5)
                        }
                        .buttonStyle(.borderedProminent)

                        Button("Reset") {
                            simulation.reset()
                        }
                        .buttonStyle(.bordered)
                        .tint(.red)

                        VStack(alignment: .leading, spacing: 4) {
                            HStack {
                                Text("Speed")
                                Spacer()
                                Text("\(simulation.speed)x")
                                    .foregroundStyle(.secondary)
                            }
                            .font(.caption2)

                            Slider(value: Binding(
                                get: { Double(simulation.speed) },
                                set: { simulation.speed = Int($0) }
                            ), in: 1...20, step: 1)
                        }
                        .frame(width: 140)
                    }
                    .padding(.horizontal, 12)
                    .padding(.vertical, 8)
                    .background(.ultraThinMaterial, in: RoundedRectangle(cornerRadius: 10))
                }
                .padding(.top, 12)
                .padding(.horizontal, 12)
            }
            .onReceive(timer) { _ in
                simulation.stepFrame()
            }
        }
    }
}

private struct StatView: View {
    let label: String
    let value: String
    let color: Color

    var body: some View {
        HStack(spacing: 4) {
            Text(label)
                .foregroundStyle(.secondary)
                .font(.caption2)
            Text(value)
                .foregroundStyle(color)
                .font(.headline)
        }
    }
}

private struct LegendItem: View {
    let color: Color
    let label: String

    var body: some View {
        HStack(spacing: 4) {
            Circle()
                .fill(color)
                .frame(width: 8, height: 8)
            Text(label)
                .font(.caption2)
        }
    }
}

private extension Color {
    static let emerald = Color(red: 0.2, green: 0.8, blue: 0.55)
    static let amber = Color(red: 0.96, green: 0.74, blue: 0.2)
}

#Preview {
    ContentView()
}
